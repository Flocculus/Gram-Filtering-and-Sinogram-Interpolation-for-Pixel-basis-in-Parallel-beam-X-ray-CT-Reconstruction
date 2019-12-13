function HTg = adjoint_space(g, yStart, yStep, ySize, Ps, ...
	xStart, xStep, xSize, phi)
% function to compute the x-ray adjoint cTilde = H'g in the space domain in
% 2D or 3D
%
% input:
% g - measurements
%
% Mike McCann 2015

numThetas = size(Ps,3);


%% loop
HTg = zeros(xSize);

xkGrid = makeGrid(xStart, xStep, xSize);
ymGrid = makeGrid(yStart, yStep, ySize);

for thetaInd = 1:numThetas
	gIndOffset = (thetaInd-1)*prod(ySize);
	
	P = Ps(:,:,thetaInd);
	offset = P * xkGrid;
		
	for yInd = 1:size(ymGrid,2)
		args = abs(bsxfun(@minus, ymGrid(:, yInd), offset));
		kernelVals = phi.xray(args, P);
		HTg(:) = HTg(:) + kernelVals' .* g(yInd + gIndOffset );
	end
	
end

global debugFlag
if debugFlag
	% plot the true values of the adjoint
	a = zeros(numel(HTg),1);
	for yInd = 1:size(ymGrid,2)
		args = abs(bsxfun(@minus, ymGrid(:, yInd), offset));
		kernelVals = phi.xray(args, P);
		a(:) = a(:) + kernelVals' .* g(yInd + gIndOffset );
	end
	figure
	hold on
	vals = sortrows([offset' a]);
	plot(vals(:,1), vals(:,2))
	
end