classdef SincFunc < Func
% phi(x) = sin( 1/2 w0 x ) / ( 1/2 w0 x) = sin(pi x) / (pi x)
%
%
%  "Foundations of Signal Processing" Example 4.5, etc.
%
% sinc is just j_0(z) = sqrt(pi/2*z) * J_{1/2}(z)
properties
	w0 = 2*pi
	
end

methods
	function obj = SincFunc()
		obj.freqRadius = 1/2 * obj.w0; % = pi
		
		obj.L1Norm = 2*pi/obj.w0;
		obj.L2Norm = 2*pi/obj.w0;
		obj.D = 1;
		
	end
	
	function y = eval(obj, x)
		y = sin( .5 * obj.w0 * x) ./ ( .5 * obj.w0 * x);
		y(x==0) = 1;
		
	end
	
	function func = abelTransformFunc(this)
		% x-ray projection of the isotropic version
		% Poularikas A. D. "The Radon and Abel Transform"
		% The Handbook of Formulas and Tables for Signal Processing.
		% Table 16.2
		% sinc = sin(pi x) / x
		% Abel transform of sinc(2*a*r) = 1/(2a) J_0 (2*pi*a*x)
		
		c = this.L1Norm;
		m = 1 / pi;
		func = IsotropicGeneratingFunc(BesselJFunc(), 1, c, m);

		
	end
	
	function p = integrate(this, y)
		% integral from -inf to y of the phi
		
		% integral over all time of sin(x)/x is pi which explains the
		% additive offset
		
		p = (2/this.w0) * (pi/2  + siNumeric( 1/2 * y * this.w0 ) );%%%%%%%%%%%%%%%%%%%%%%siNumeric

		
	end
	
	function val = autocorrelation(this, x)
		val = this.eval(x);
	end
	
	% fourier transformph
	% "Foundations of Signal Processing" Table 4.1
	
	
end

end