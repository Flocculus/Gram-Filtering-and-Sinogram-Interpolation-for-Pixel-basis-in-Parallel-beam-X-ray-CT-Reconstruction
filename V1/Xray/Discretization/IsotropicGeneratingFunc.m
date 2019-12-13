classdef IsotropicGeneratingFunc < GeneratingFunc
% ND isotropic generating function

properties
	phi1D
end

methods
	function obj = IsotropicGeneratingFunc(phi1D, D, c, m)
		obj.D = D;
		obj.phi1D = phi1D;
			
		if ~exist('c', 'var') || isempty(c)
			obj.c = 1;
		elseif length(c) ~= 1
			error('c must be a scalar');
		else
			obj.c = c;
		end
		
		if ~exist('m', 'var') || isempty(m)
			obj.m = 1;
		elseif length(m) ~= 1
			error('m must be a scalar');
		else
			obj.m = m;
		end
		
		obj.spaceRadius(1:D) = obj.phi1D.spaceRadius * m;
		obj.freqRadius(1:D) = obj.phi1D.freqRadius / m; % todo: fix
		
	end
	
	function val = eval(obj, x)
		
		val = sqrt( sum( x.^2, 1 ) );
		
		val = obj.c * obj.phi1D.eval(val / obj.m); 
	end

	
	function val = integrate(obj, x)
		% integral from -infinity to x_d in each dimension, can call on a 
		% D x N matrix, the result is 1 X N.
		if obj.D ~=1
			error('multi-D integrals not done yet')
		end
		
		val = obj.c * obj.m * obj.phi1D.integrate(x/obj.m);
	end
	
	function val = autocorrelation(obj, x)
		xNorm = sqrt( sum( x.^2, 1 ) );
		val = obj.c^2 * obj.m * obj.phi1D.autocorrelation(xNorm/obj.m);
	end
	
	
end

end