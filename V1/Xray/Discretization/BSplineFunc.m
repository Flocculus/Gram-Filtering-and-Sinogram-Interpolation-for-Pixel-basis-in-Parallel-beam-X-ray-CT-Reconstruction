classdef BSplineFunc < Func
	properties (SetAccess = immutable)
		degree % degree
	end
	
	properties (Access = private)
		coeffs
		
	end
	
	
	methods
		function obj = BSplineFunc(n)
			obj.degree = n;
			obj.spaceRadius = (n + 1) / 2;
			obj.D = 1; % # of dimensions
			
		end
		
		function y = eval(obj, x)
			% from "splines a perfect fit", box 1 (10), but computed in a
			% smart way. We track the coefficients of the piecewise polynomial as k
			% changes. The benefit of this (over naive eq 10) is that the B-spline is exactly equal to
			% zero outside its support.
			
			if obj.degree == 0
				y = double(abs(x) < .5);
				y(abs(x) == 1/2) = .5;
				return
			end
			
			y = zeros(size(x));
			
			if isempty(obj.coeffs)
				obj.coeffs = zeros(obj.degree+1, obj.degree+1);
				
				for k = 0 : obj.degree
					
					newCoeffs = [-k+(obj.degree+1)/2 1];
					for i = 2:obj.degree
						newCoeffs = conv(newCoeffs,  [-k+(obj.degree+1)/2 1]);
					end
					
					if k >= 1
						obj.coeffs(k+1, :) = obj.coeffs(k, :) + nchoosek(obj.degree+1,k) * (-1)^k * newCoeffs;
					else
						obj.coeffs(k+1, :) =  nchoosek(obj.degree+1,k) * (-1)^k * newCoeffs;
					end
				end
				
			end
			
			
			for k = 0:obj.degree
				l = k-(obj.degree+1)/2;
				
				for i = 1:size(obj.coeffs, 2)
					y( x >= l & x < l+1 ) = y( x >= l & x < l+1 ) ...
						+ x( x >= l & x < l+1 ).^(i-1) * obj.coeffs(k+1, i);
				end
			end
			y =  1/factorial(obj.degree) * y;
			
		end
		
		function y = integrate(this, x)
			% integral from -inf to y of the phi
			
			Bint = BSplineFunc(this.degree+1);
			
			kMin = max( ceil((min(x) - .5) - Bint.spaceRadius), 0);
			kMax = floor( (max(x) - .5) + Bint.spaceRadius );
			
			% thus x - 1/2 - k > -Bint.spaceRadius
			
			y = zeros(size(x));
			for k = kMin:kMax
				y = y + Bint.eval( x - .5 - k);
			end
			
			
		end
		
	function func = abelTransformFunc(this)
		if this.degree ~= 0
			error('abel transform only implemented for degree 0');
		end
		
		func = @(x) 2 * sqrt( .25 - x.^2 ) .* double(abs(x) <= .5);
		func = IsotropicGeneratingFunc(GenericFunc(func), 1, 1, 1);
		
	end
		
	end
end