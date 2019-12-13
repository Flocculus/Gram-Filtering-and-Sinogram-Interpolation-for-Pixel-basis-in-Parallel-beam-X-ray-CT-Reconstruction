classdef BSplineKernel < SeparableGeneratingFunc
	
	properties (SetAccess = immutable)
		degree % degree
	end
	
	methods
		function this = BSplineKernel(n, D, c, m)
			
			[funcs{1:D}] = deal(BSplineFunc(n));
			
			this = this@SeparableGeneratingFunc(funcs, c, m);
			
			this.degree = n;
		end
	
		
		
		function y = reconstruct(obj, c, x, boundaryCond)
			if ~exist('boundaryCond', 'var') || isempty(boundaryCond)
				boundaryCond = 'zero';
				% todo: handle boundary conditions here?
			end
			
			D = size(x,1);
			
			if D == 1;
				c = c(:);
			end

			switch boundaryCond
				case 'mirror'
					% (default in interpol)
					if D == 1
						y = interpol(c, x, obj.degree);
					elseif D == 2
						y = interpol(c, x(1,:), x(2,:), obj.degree);
					end
				case 'zero'
					% only these x vals are nonzero
					inside = all(x >= -obj.spaceRadius) ...
					& all(bsxfun(@le, x, size(c)' - 1 + obj.spaceRadius));
					y = zeros(1, size(x,2));
					
					% have to pad c to circumvent mirror boundaries in
					% interpol
					pad =  floor(obj.spaceRadius);
					
					if D == 1
						y(inside) = interpol(padarray(c,pad), x(inside) +pad , obj.degree);
					elseif D == 2
						y(inside) = interpol(padarray(c,[pad pad]), x(1,inside) +pad , x(2,inside)+pad, obj.degree);
					end
			end
		end
		
		function  y = density_estimate(obj, N, x)
			c = zeros(N,1);
			c(1) = 1;

			
			% only these y vals are nonzero
			inside = x >= -obj.spaceRadius & x <= length(c)-1+obj.spaceRadius;
			
			% have to pad c to circumvent mirror boundaries in
			% interpol
			pad =  floor(obj.spaceRadius);
			y = density_estimate(padarray(c,2*pad, 'post'), x(inside) , obj.degree);
			y = y(1:N);
		end
		
		
	end % methods
		
		
end % classdef