classdef SeparableGeneratingFunc < GeneratingFunc
% phi(x) = prod_{d=0}^{D-1} c_d phi_d( x_d / m_d ) 
%
properties
	phis
end

methods
	function obj = SeparableGeneratingFunc(phis, c, m)
		% cell array of Funcs
		
		Ds = cellfun(@(x) x.D, phis);
		if any(Ds ~= 1)
			error('assuming 1D component functions');
		end
		
		obj.D = length(phis);
		
		obj.L1Norm = prod(m);
		
		if ~exist('c', 'var') || isempty(c)
			obj.c = ones(1,D);	
		elseif length(c) == 1
			obj.c(1:obj.D) = c;
		elseif length(c) == obj.D
			obj.c = c;
		end
		
		if ~exist('m', 'var') || isempty(m)
			obj.m = ones(1,D);
		elseif length(m) == 1
			obj.m(1:obj.D) = m;
		elseif length(m) == obj.D
			obj.m = m;
		end
		
		
		obj.phis = phis;
		for i = 1:obj.D
			obj.spaceRadius(i) =  obj.phis{i}.spaceRadius * obj.m(i);
			obj.freqRadius(i)  =  obj.phis{i}.freqRadius  / obj.m(i);
		end
		
		% the x-ray transform of an isotropic function is also isotropic,
		% and its profile is the Abel transform of the original profile
		
	end
	
	function val = eval(obj, x)
		obj.checkInput(x);
		
		val = ones(1, size(x, 2));
		for d = 1:obj.D
			val = val * obj.c(d) .* (obj.phis{d}.eval( x(d,:) / obj.m(d) ));
		end
		val = val;
	end
	
	
	function val = integrate(this, x)
		% integral from -infinity to x_d in each dimension, can call on a
		% D x N matrix, the result is 1 X N.
		val = ones(1, size(x, 2));
		for d = 1:this.D
			val = val .* this.integrate1D(x(d,:), d);
		end
		
	end
	
		function val = integrate1D(this, x, d)
		% integral from -infinity to x_d in each dimension, can call on a
		% D x N matrix, the result is 1 X N.

		val = this.m(d) * this.c(d) .* this.phis{d}.integrate(x/this.m(d));

		
	end
	
	
	function val = autocorrelation(this, x)
		val = ones(1, size(x, 2));
		for d = 1:this.D
			val = val * this.c(d) .* (this.phis{d}.autocorrelation( x(d,:) / this.m(d) ));
		end
	end
	
end

methods (Static)
	function test(phi)
		dx = .1;
		
		if all(isfinite(phi.spaceRadius))
			x = -phi.spaceRadius(1):dx:phi.spaceRadius(1);
			y = -phi.spaceRadius(2):dx:phi.spaceRadius(2);
		else
			x = -20:dx:20;
			y = x;
		end
		
		[X, Y] = ndgrid(x, y);
		Z = reshape( phi.eval([X(:)'; Y(:)']), size(X) );
		figure
		imagesc(Z, 'xdata', x, 'ydata', y);
		figure
		h = surf(X, Y, Z);
		
		
		set(h, 'edgecolor','none')
	end
	
end

end