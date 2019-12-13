% method to define all the enums that are used throughout

METHODS = enum({'exact', 'standard', 'orthogonal', 'oblique','boxspline'});
PHANTOMS = enum({'ellipse', 'tinyElls', 'bigSinc', 'box', 'bigKBW',...
	'smallKBW', 'multiKBW', 'noise', 'smallSinc', 'sheppSmooth', 'filament', 'SheppLogan', 'forbuild', "realImg"});