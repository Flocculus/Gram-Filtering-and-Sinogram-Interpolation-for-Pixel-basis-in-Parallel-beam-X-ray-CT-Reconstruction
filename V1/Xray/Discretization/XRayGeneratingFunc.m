classdef (Abstract) XRayGeneratingFunc < handle
	% 2D/3D generating function for X-ray reconstruction
	
	methods (Abstract)
		val = xray(obj, y, P);
		% P corresponds to P_{\theta^\perp} in "Fast 3D Reconstruction Method
		% for Differential Phase Contrast X-ray", no transpose. So the
		% basis is in the rows
		
	end
	
	
	methods (Static)
		function testXray(phi)
			dx = .05;
			
			t = pi/6;
			%t = 0;
			
			if all(isfinite(phi.spaceRadius))
				x = -phi.spaceRadius(1):dx:phi.spaceRadius(1);
				y = -phi.spaceRadius(2):dx:phi.spaceRadius(2);
			else
				x = -100:dx:100;
				y = x;
			end
			
			[X, Y] = ndgrid(x, y);
			Z = reshape( phi.eval([X(:)'; Y(:)']), size(X) );
			
			
			Z = imrotate(Z, -t/pi*180, 'crop');
			figure
			hold on
			plot(x, sum(Z,1) * dx)
			

			plot(x, phi.xray(x, [cos(t) sin(t)]))
			legend({'numerical', 'analytical'})
			
		end
		
		
	end
	
end


	