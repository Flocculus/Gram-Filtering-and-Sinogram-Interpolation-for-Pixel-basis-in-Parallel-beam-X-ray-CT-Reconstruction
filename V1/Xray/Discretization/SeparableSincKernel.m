classdef SeparableSincKernel < XRayGeneratingFunc & SeparableGeneratingFunc
	% 2D or 3D separable sinc kernel.
	% for 3D, only fixed axis projections are allowed.
	properties
		
		xrayTransform;
	end
	
	methods
		function this = SeparableSincKernel(D, m, c)
			
			if length(m) == 1;
				m = ones(1,D)*m;
			elseif length(m) == D
			else
				error('wrong m length');
			end
			
			if ~exist('c', 'var') || isempty(c)
				c = ones(1,D);
			end
			
			[funcs{1:D}] = deal(SincFunc());
			this = this@SeparableGeneratingFunc(funcs, c, m);

			if this.D == 2
				this.xrayTransform = IsotropicGeneratingFunc(SincFunc(), this.D-1, 1, 1);
			elseif this.D == 3
				this.xrayTransform = SeparableGeneratingFunc({SincFunc() SincFunc()}, [1 this.c(3)], [1 this.m(3)]);
			end
		end
		
		
		function val = xray(this, y, P)
			% from "Intro to sparse stochastic processes: 10.3.4
			
			xrayPhi = this.xrayTransformFunc(P);
			val = xrayPhi.eval(y);
			
		end
		
		function func = xrayTransformFunc(this, P)
			% if given empty P, assume the projection is along the x_0 axis
			if isempty(P)
				P = [zeros(this.D-1, 1) eye(this.D-1) ];
			end

			
			if this.D == 2
				m = max(abs(P) .* this.m);
				c = prod(this.c) * this.L1Norm / m;
				
			elseif this.D == 3
				m = max(abs(P(1,1:2)) .* this.m(1:2)); %assume fixed axis
				c = prod(this.c(1:2)) * this.L1Norm / m;	
			end
			this.xrayTransform.c(1) = c;
			this.xrayTransform.m(1) = m;
			
			func = this.xrayTransform;
			
			
		end
		
	end
	
	
end