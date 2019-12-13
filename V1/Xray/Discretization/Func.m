classdef (Abstract) Func < handle
	
properties 
	D; % number of dimensions
	
	L1Norm;
	L2Norm;
	
	spaceRadius = inf;
	freqRadius = inf;
end

methods(Abstract)
	val = eval(x)
end

methods
	function val = hat(w)
		error('no Fourier transform implemented')
	end
	
	% would be cool to overload () for Func's, but can't get it to work
	% in cases like phi.phi1D(0)
	% 	function varargout = subsref(obj, S)
	% 		switch S(1).type
	% % 			case '()'
	% % 				varargout{1} = obj.eval(S.subs{1});
	% 			otherwise
	% 				[varargout{1:nargout}] = builtin('subsref',obj,S);
	% 		end
	% 	end
	
	function val = autocorrelation(this, ~)
		val = NaN;
	end
	
	function func = abelTransformFunc(this)
		% useful because the abel transform of a 1D function is exactly the
		% x-ray transform of the isotropic kernel formed from it
		error('no abel transform implemented');
	end
	
end

methods (Static)
	function test_integral(f)
		R = min(f.spaceRadius,20);
		dx = .5;
		x = -R:dx:R;
		[integrals{1:2}] = deal(zeros(size(x)));
		% numerical
		for i = 1:length(x)
			integrals{1}(i) = integral(@(x) f.eval(x), -10000, x(i));
		end
		% analytical
		integrals{2} = f.integrate(x);
		
		figure
		plot(x, integrals{1}, x, integrals{2})
		
	end
end
	
end