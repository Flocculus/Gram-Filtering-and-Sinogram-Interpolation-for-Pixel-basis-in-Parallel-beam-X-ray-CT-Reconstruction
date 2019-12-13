classdef KBWFunc < Func
% 1D generalized kaiser bessel window
properties
	m % bessel order
	alpha % smoothness
	a % radius of nonzero support
end

methods
	function obj = KBWFunc(m, alpha, a)
		
		obj.m = m;
		obj.alpha = alpha;
		obj.a = a;
		obj.D = 1;
		
		obj.spaceRadius = obj.a;
		
		% based on the FT at zero
		A = (2*pi)^(obj.D/2) * obj.a^obj.D * obj.alpha^obj.m / besseli(obj.m, obj.alpha);
		obj.L1Norm = A * besseli( obj.D/2+obj.m, obj.alpha ) / obj.alpha.^(obj.D/2+obj.m);
		
	end
	
	function p = eval(obj, x)
		% x - D X N matrix of inputs
		%
		% "Multidimensional digital image representations using generalized
		% Kaiser-Bessel window functions" (A1)
		

		xSq = x.^2;

		p = zeros(size(xSq));
		%
        %if xSq<=obj.a ^2
        %    isNonezero = 1;
        %else
        %    isNonezero = 0;
        %end
        %
		isNonzero = xSq <= obj.a ^2;
		A = sqrt( 1 - xSq(isNonzero) / obj.a^2 );
		p(isNonzero) = A.^obj.m .* besseli(obj.m, obj.alpha*A) / besseli(obj.m, obj.alpha);
		
	end
	
	function func = abelTransformFunc(this)
		% x-ray projection of generalized kaiser bessel window
		% "Multidimensional digital image representations using generalized
		% Kaiser-Bessel window functions" (A7)
		
		c =  this.a / besseli(this.m, this.alpha) * ...
			sqrt(2*pi / this.alpha) * besseli(this.m + .5, this.alpha);
		
		func = IsotropicGeneratingFunc(KBWFunc(this.m+.5, this.alpha, this.a), 1, c, 1);
		
	end
	
	
	function g = grad(obj, x, dim)
		% gradient of the kbw at the points specified in x along dimension
		% dim.
		% s - obj.D X number of points to evaluate
		% dim - scalar, {1, 2, 3}, which dimension to differentiate along
		%
		% "Multidimensional digital image representations using generalized
		% Kaiser-Bessel window functions" (A11)
		xSumSq = sum( x.^2, 1);
		
		isNonzero = xSumSq <= obj.a ^2;
		
		z = obj.alpha * sqrt( 1 - xSumSq / obj.a^2 );
		
		A = -(1 / obj.a) / (obj.alpha^(obj.m-2) * besseli(obj.m,obj.alpha));
		
		c1 = zeros(size(xSumSq));
		c1( isNonzero) = x(dim, isNonzero) / obj.a .* z(isNonzero).^(obj.m-1) .*besseli(obj.m-1, z(isNonzero));
		
		g = A * c1;
	end
	
	
	
	
	function wHat = hat(obj, r)
		% NORMALIZED fourier transform of generalized kaiser bessel window
		
		% "Multidimensional digital image representations using generalized
		% Kaiser-Bessel window functions" (A3)
		
		s = r;
		
		
		A = (2*pi)^(obj.D/2) * obj.a^obj.D * obj.alpha^obj.m / besseli(obj.m, obj.alpha);
		
		wHat = zeros(size(s));
		
		c1 = abs(2*pi*obj.a*s) < obj.alpha;
		B =  sqrt( obj.alpha^2 - (2*pi*obj.a*s(c1)).^2 );
		wHat(c1) = A * besseli( obj.D/2+obj.m, B ) ./ B.^(obj.D/2+obj.m);
		
		
		c2 = abs(2*pi*obj.a*s) >= obj.alpha;
		C =  sqrt( (2*pi*obj.a*s(c2)).^2 - obj.alpha^2);
		wHat(c2) =  A * besselj( obj.D/2+obj.m, C ) ./ C.^(obj.D/2+obj.m);
		
	end
	
	function wHat = xrayHat(obj, r)
		% NORMALIZED fourier transform of projection of generalized kaiser bessel window
		% that is, r=1 -> w = 2*pi
		% "Multidimensional digital image representations using generalized
		% Kaiser-Bessel window functions" (A3)
		
		s = sqrt( sum( r.^2, 1) );
		
		m1 = obj.m + 0.5;
		n1 = obj.D - 1;
		
		A = (2*pi)^(n1/2) * obj.a^n1 * obj.alpha^m1 / besseli(m1, obj.alpha);
		
		wHat = zeros(size(s));
		
		c1 = abs(2*pi*obj.a*s) < obj.alpha;
		B =  sqrt( obj.alpha^2 - (2*pi*obj.a*s(c1)).^2 );
		wHat(c1) = A * besseli( n1/2+m1, B ) ./ B.^(n1/2+m1);
		
		
		c2 = abs(2*pi*obj.a*s) >= obj.alpha;
		C =  sqrt( (2*pi*obj.a*s(c2)).^2 - obj.alpha^2);
		wHat(c2) =  A * besselj( n1/2+m1, C ) ./ C.^(n1/2+m1);
		
		wHat = wHat *  besseli(m1, obj.alpha) *  obj.a ./ besseli(obj.m, obj.alpha) .* sqrt(2*pi/obj.alpha);
	end
	
	
end

methods(Access = protected)

	
	
end



methods(Static)
	function test()
		r = 1;
		alpha = 10.83; % 7.91 or 10.83 are optimal
		m = 2;
		
		phi = KBWFunc(m, alpha, r, 2);
		
		x = linspace(-r,r, 101);
		xStep = x(2)-x(1);
		[x0, x1] = ndgrid(x);
		
		phiVals = reshape(phi.eval([x0(:) x1(:)]'), size(x0));
		imagesc(x, x, phiVals')
		colorbar
		axis xy
		xlabel('x_0');
		ylabel('x_1');
		
		
		phiSum = sum(phiVals(:)) * (xStep)^2
		phiHatAt0 = phi.hat(0)
		
		phiProjVals = phi.xray(x);
		phiProjSum = sum(phiProjVals) * xStep
		phiProjHatAt0 = phi.xrayHat(0)
	end
	
	
	function vals = kbwLookup(x, table, width)
		% unused code, might be useful later
		vals = zeros(size(x,1), size(x,2));
		inds = round(2*abs(x)/width * size(table,1)) + 1;
		isGood = inds <= size(table,1);
		vals(isGood) = table(inds(isGood));
	end
end
end