classdef XRay < LinOp
	properties (SetAccess = private, GetAccess = public)
		D
		
		xStart
		xStep
		xSize
		
		yStart
		yStep
		ySize
		Ps
		
		phi
		numProj
		
		yUpsampleRateHtH = 10;
		
		HTHHat
		HTHkernel
		
		norm1Kernel
		
		options
	
		
	end
	
	methods
		function this = XRay(xStart, xStep, xSize, yStart, yStep, ySize, Ps, phi, varargin)
			this.D = size(xStart,2);
			this.phi = phi;
			this.Ps = Ps;
			this.numProj = size(Ps,3);
			
			
			this.xStart = xStart;
			this.xStep = xStep;
			this.xSize = xSize;
			this.sizein = xSize;
					
			
			this.yStart = yStart;
			this.yStep = yStep;
			this.ySize = ySize;
			this.sizeout = [ySize this.numProj];
			
			% parse other options
			p = inputParser;
			addOptional(p,'adjointMethod', 2); % LUT
			addOptional(p,'adjointDegree', 1); % 
			addOptional(p, 'adjointUpsampleRate', 2);
			
			addOptional(p, 'kernelMethod', 'analytical');
			parse(p,varargin{:});
			this.options = p.Results;
			
			%% check consistency of options
			
			% kernel method
			if strcmp(this.options.kernelMethod, 'analytical')
				xrayPhi = this.phi.xrayTransformFunc([]);
				if isnan(xrayPhi.autocorrelation(yStart'))
					warning('analytical kernel method changed to discreteConv because xrayPhi autocorrelation is unknown');
					this.options.kernelMethod = 'discreteConv';
				end
			end
		end
		
		function Hc = Apply(this, c) % Apply the operator
			if isa(this.phi, 'SeparableSincKernel')
				assert(all(this.phi.c==1), 'we do not handle scaled sinc kernels');
				Hc = sinc_grid_forward(this.xStart, this.xStep, this.xSize, ...
					this.yStart, this.yStep, this.ySize, this.Ps, this.phi.m, c);
			else
				Hc = forward_project(c, this.yStart, this.yStep, this.ySize, this.Ps, ...
					this.xStart, this.xStep, this.xSize, this.phi);
			end
		end
		
		function c = FBP(this, g)
			padSize = 3*size(g,2);
			gPad = padarray(g, [padSize 0]);
			
			gHat = fft(gPad, [], 1);
			om = (0:size(gHat,1)-1)/size(gHat,1) * 2 * pi;
			om(om>pi) = 2*pi - om(om>pi);
			om = om / pi; 
			
			gFiltHat = bsxfun(@times, gHat, om');
			gFilt = real(ifft(gFiltHat, [], 1));
			
% 			gFilt = bsxfun(@minus, gFilt, mean(gFilt([padSize+1 end-padSize],:),1) );
			
			
 			c = this.Adjoint(gFilt(padSize+1:end-padSize, :));
			c = c / prod(this.xStep) / size(g,2) * pi /2; 
		end
		
		function HTg = Adjoint(this,g) % Apply the adjoint
			degree = this.options.adjointDegree;
			upsampleRate = this.options.adjointUpsampleRate;
			
			args = {g, this.yStart, this.yStep, this.ySize, this.Ps, this.xStart, this.xStep, this.xSize, this.phi};
			switch this.options.adjointMethod
				case 1
					HTg = adjoint_space(args{:});
				case 2
					HTg = adjoint_LUT(args{:}, upsampleRate, degree);
				case 3
					HTg = adjoint_oblique(args{:}, upsampleRate, degree, false);
				case 4
					HTg = adjoint_oblique(args{:}, upsampleRate, degree, true);
				case 5
					HTg = adjoint_LUT_noSmooth(args{:}, upsampleRate, degree);
				case 6
					HTg = adjoint_oblique_noSmooth(args{:}, upsampleRate, degree, false);
				case 7
					HTg = adjoint_oblique_noSmooth(args{:}, upsampleRate, degree, true);
			end
			
		end
		
		function precomputeHtH(this)
			if isempty(this.HTHkernel)
	
				
				
				kSize = (this.xSize-1)*2+1; % ensure odd so we can sample zero
				kStep = this.xStep;
				kStart = -(kSize-1)/2 .* kStep;
				kGrid = makeGrid(kStart, kStep, kSize);
				
				switch this.options.kernelMethod
					case 'analytical'
						% compute the autocorrelation of the x-ray
						% transform of the basis function analytically,
						% which requires that you know it.
						h = zeros([kSize 1]);
						for projInd = 1:this.numProj
							if this.D == 3
								fprintf('%d, ', projInd);
							end
							P = this.Ps(:,:,projInd);
							xrayPhi = this.phi.xrayTransformFunc(P);
							h(:) = h(:) + xrayPhi.autocorrelation( P * kGrid )';
						end
						if this.D == 3
							fprintf('\n')
						end
						
						h = 1/prod(this.yStep) * h; % comes from the math!!
						
					case 'discreteConv'
						% compute the autocorrelation of the x-ray
						% transform of the basis function using an
						% upsampled discrete convolution
						yUpStep =  ones(1, this.D-1)/2; % 2x upsampling
						doMex = this.D==3;
						isIsotropic = false;
						
						
						h = HTH_kernel_LUT(yUpStep, this.Ps, ...
							kStart, kStep, kSize, this.phi, doMex, isIsotropic);
						
					case 'forwardBackward'
						% calculate the sinogram of the impulse, then apply
						% the adjoint to it. This is exact if the exact
						% adjoint is used
						yGrid = makeGrid(this.yStart, this.yStep, this.ySize);
						
						gInds = make_grid_vectors(ones(1,this.D-1), ones(1,this.D-1), this.ySize);
						g = zeros([this.ySize size(this.Ps,3)]);
						for projInd = 1:this.numProj
							g(gInds{:}, projInd) = reshape( this.phi.xray(yGrid, this.Ps(:,:,projInd)), [this.ySize 1] );
						end
						
						degree = this.options.adjointDegree;
						upsampleRate = this.options.adjointUpsampleRate;
						
						args = {g, this.yStart, this.yStep, this.ySize, this.Ps, kStart, kStep, kSize, this.phi};
						switch 2
							case 1
								h = adjoint_space(args{:});
							case 2
								h = adjoint_LUT(args{:}, upsampleRate, degree);
							case 3
								h = adjoint_oblique(args{:}, upsampleRate, degree, false);
							case 4
								h = adjoint_oblique(args{:}, upsampleRate, degree, true);
							case 5
								h = adjoint_LUT_noSmooth(args{:}, upsampleRate, degree);
							case 6
								h = adjoint_oblique_noSmooth(args{:}, upsampleRate, degree, true);
						end
						
						% for some kernels, could use the fact that it is
						% symmetric and only calculate part. or force it to
						% be symmetric with this code;
						% h = (h + rot90(h,2))/2;
						
				end
				
				
				this.HTHkernel = h;
				this.norm1Kernel = sum(h(:));
				% norm1Kernel is approximately how much a delta in c adds to HTg
				
				% handle the padding of hIn
				this.HTHHat = fftn( circshift(h, -this.xSize +1 ) );
				%this.norm1Kernel = this.HTHHat(1);
			end
		end
		
		function g = HtH(this, c)
			
			this.precomputeHtH() % compute the kernel if needed

			g = convxh(c, this.HTHHat, [], true); 
			
			%g = this.Adjoint( this.Apply( c ));
			
			
		end
		
		function y = HHt(this,x) %  Apply the HHt matrix
			if this.issquare   % HtH =HHt
				y = this.HtH(x);
			else
				y = this.Apply(this.Ajoint(x));
			end
		end
		
		function y = HtWH(this,x,W) %  Apply the HtH matrix
			if (isscalar(W) && isreal(W))
				y = this.HtH(x);
			else
				assert(isa(W,'LinOp'),'W must be a LinOp');
				y = this.Adjoint(W.Apply(this.Apply(x)));
			end
		end
		
		function Inverse(this,~) % Apply the inverse
			if this.isinvertible
				error('Inverse not implemented');
			else
				error('Operator not invertible');
			end
		end
		
		function AdjointInverse(this,~) % Apply the inverse
			if this.isinvertible
				error('AdjointInverse not implemented');
			else
				error('Operator not invertible');
			end
		end
	end

		
	
	
	
end
