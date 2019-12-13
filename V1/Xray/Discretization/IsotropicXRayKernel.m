classdef IsotropicXRayKernel < IsotropicGeneratingFunc & XRayGeneratingFunc
	properties
		xrayPhi
	end
	
	
	
	methods
		function this = IsotropicXRayKernel(phi1D, D, c, m)
			if D < 2
				error('you are asking for 1D X-ray tomography');
			end
			if ~exist('c', 'var') || isempty(c)
				c = 1;
			end
			if ~exist('m', 'var') || isempty(m)
				m = 1;
			end
			
			
			this = this@IsotropicGeneratingFunc(phi1D, D, c, m);
			
			% the x-ray transform of an isotropic function is also isotropic,
			% and its profile is the Abel transform of the original profile
			

			xrayPhi1D = phi1D.abelTransformFunc();
			xrayPhi1D.c = this.c * xrayPhi1D.c * this.m;
			xrayPhi1D.m = this.m * xrayPhi1D.m;
			this.xrayPhi = IsotropicGeneratingFunc(xrayPhi1D, D-1, 1, 1);
		end
		
		function val = xray(this, y, ~)
			val = this.xrayPhi.eval( y );
			
		end
		
		function func = xrayTransformFunc(this, ~)
			func = this.xrayPhi;
			
		end
		
	end
	
	methods (Static)
	
		
	end
	
end