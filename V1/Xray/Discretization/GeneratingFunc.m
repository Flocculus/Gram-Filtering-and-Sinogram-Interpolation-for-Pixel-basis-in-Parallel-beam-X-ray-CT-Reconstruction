classdef GeneratingFunc < Func
% ND generating function: c phi( x / m )

properties	
	c; % multiplicative constant
	m; % scaling
end

methods

	function checkInput(obj, x)
		if size(x,1) == obj.D
			return;
		else
			error('wrong input size, should be %d by N', obj.D);
		end
	end
	
	function y = conv(this, xStart, xStep, xSize, x, yStart, upRate, ySize)
		% function y = conv(this, xStart, xStep, xSize, x, yStart, yUpsampleRate, ySize)
		%
		% compute the convolution between phi and the values in x, return
		% results at the locations specified in the y grid
		
		if length(upRate) == 1
			upRate(1:this.D) = upRate;
		end
		
		
		% the k grid is where the convolution happens
		kStep = xStep ./ upRate;
		
		xEnd = xStart + xStep .* (xSize-1);
		yEnd = yStart + kStep .* (ySize-1);
		
		left = min([xStart; yStart]);
		right = max( [xEnd; yEnd] );
		
		
		kStart = xStart + floor((left-xStart)./kStep) .* kStep;
		kSize = ceil((right - kStart)./kStep) + 1;

		% upsample x
		xUp = zeros( [kSize 1] );
		upVects = make_grid_vectors(round((xStart - kStart)./kStep)+1, upRate, xSize);
		xUp(upVects{:}) = x;
		
		% get phi values on the same grid
		yOff = ((yStart - xStart) ./ kStep - floor( (yStart - xStart) ./ kStep )) .* kStep;
		r = min( [floor(this.spaceRadius ./ kStep); kSize-1]);
		rStart = -r;
		rSize = 2*r + 1;
		
		rGrid = makeGrid(rStart, ones(size(rSize)), rSize);
		phiVals = reshape( this.eval( bsxfun(@plus, bsxfun(@times, rGrid,  kStep'), yOff') ), [rSize 1]);
		
		
		% do the convolution
		yAll = convxh(xUp, phiVals, rSize); 
		
		yVectors = make_grid_vectors( round((yStart-yOff-kStart)./kStep) + 1, ones(size(ySize)), ySize);
		y = yAll(yVectors{:});
		
		
	end
	
	
	
	function f = reconstruct(obj, xStep, c)
		xSize = size(c);
		extraDims = prod(xSize(obj.D+1:end));
		xSize = xSize(1:obj.D);
		
		if extraDims > 1
			error('wrong c dimension');
		end
		if length(xStep) ~= obj.D
			error('wrong xStep length');
		end
		
		xStart = zeros(size(xStep));
		
		
		f = obj.conv(xStart, xStep, xSize, c, xStart, ones(size(xSize)), xSize);

	end
	
	
	function S = saveobj(obj)
		S = struct();
		for field = fields(obj)'
			f = field{1};
			S.(f) = obj.(f);
		end
		S.class = class(obj);
	end
	
	
	
end

methods (Static) %% tests -------------------------------------------------
	function test_getProjSeparability()
		phi = BlankTest();
		phi.D = 5;
		phi.separableGroups = [1 2 1 3 3];
		P =[...
			1 0 1 0 0;
			0 0 1 0 0;
			0 0 0 0 1;
			0 1 0 1 0];
		assert(all(phi.getProjSeparability(P) == [1 1 2 2]));
		
		P =[...
			1 1 0 0 0;
			0 0 1 1 0;
			0 0 0 0 1;
			0 1 0 1 0];
		assert(all(phi.getProjSeparability(P) == [1 1 1 1]));
		
		
		phi.separableGroups = [4 5 1 3 2];
		P =[...
			1 0 1 0 0;
			0 0 1 1 0;
			0 0 0 1 0;
			0 1 0 0 0];
		assert(all(phi.getProjSeparability(P) == [1 1 1 2]));
	end
	
	function test_getProjIsotropy()
		phi = BlankTest();
		phi.D = 3;
		phi.isotropicGroups = [1 1 2];
		
		% Ps(:,:,i) is the projection matrix from x-space to y-space (2x3)
		thetas = linspace(0, pi, 10);
		Ps = zeros(2, 3, length(thetas));
		Ps(2, 3, :) = 1;
		for tInd = 1:length(thetas)
			Ps(1, 1, tInd) = cos(thetas(tInd));
			Ps(1, 2, tInd) = sin(thetas(tInd));
		end
		
		phi.getProjIsotropy(Ps);
	end
	
	function obj = loadobj(S)
		const = str2func(S.class);
		obj = const();

		S = rmfield(S, 'class');
		for field = fields(S)'
			f = field{1};
			obj.(f) = S.(f);
		end
		
		
	end
	
end
end

