function [grid] = makeGrid(xStart, xStep, xSize)

% utility function for explicitly listing the points in the grid specified
% by (start,step,size). Useful with the KaiserBesselWindow object's
% methods and elsewhere.
%
% the output is a D X prod(xSize) matrix, which is why we don't use ndgrid
% directly

D = length(xStart);

vects = make_grid_vectors(xStart, xStep, xSize);

try
	grid = ndgridMatrix( vects{:} );
	
catch ME % no-mex way, takes >2X the memory of the eventual output
	warning('ndgridMatrix.c not compiled, defaulted to (slow) pure matlab. Try running compile.m');
	grid = cell(D,1);
	[grid{:}] = ndgrid( vects{:} );
	
	grid = cellfun(@(a)(a(:)), grid, 'uniformoutput', false)';
	grid = [grid{:}]';
end