function HTg = adjoint_LUT(g, yStart, yStep, ySize, Ps, ...
	xStart, xStep, xSize, phi, yUpsampleRate, degree)
% Compute the x-ray adjoint, HTg, in the space domain in
% 2D or 3D via lookup table as described in 

% M. T. McCann, M. Nilchian, M. Stampanon, and M. Unser,
% "Fast 3D Reconstruction Method for Differential Phase Contrast X-ray CT."
% 2016.
%
%
% input:
% g - sinogram, ySize X numProjections
% yStart - coordinates (in y-space) of g(1, i) (2D) or g(1, 1, i) (3D), 1 X D-1
% yStep - sampling step in y-space, 1 X D-1
% ySize - size(g), 1 X D-1
% Ps - projections matrix between y-space and x-space, D-1 X D X
%     numProjections
% xStart - coordinates (in x-space) of HTg(1, 1) (2D) or HTg(1, 1, 1) (3D), 1 X D
% xStep - sampling step in x-space, 1 X D
% xSize - size(HTg), 1 X D
% phi - object of type GeneratingFunc that defines the desired discretization
%     kernel for the problem
% yUpsampleRate - how much to upsample when computing the LUT. Higher is
%     more accurate but slower. 2 is a good compromise.

% Mike McCann
% michael.mccann@epfl.ch
% March, 2016
% 


%% handle input
% number of dimensions
D = length(xStart);

if ~( D == 2 || D == 3)
  error('Ps suggests this is a %d-D problem, which we do not handle', D);
end

if length(yUpsampleRate) == 1
	yUpsampleRate = ones(1,D-1) * yUpsampleRate;
end

numProj = size(Ps,3);

%% setup for main loop
% upsampling setup
yUpStep = yStep ./ yUpsampleRate;
yUpSize = ySize .* yUpsampleRate;
yUpStart = yStart;

yUpGrid = makeGrid(yUpStart, yUpStep, yUpSize);

upSampleVects = make_grid_vectors(ones(1,D-1), yUpsampleRate, ySize);
yInds = make_grid_vectors(ones(1,D-1), ones(1,D-1), ySize);

% grid to sample phi on
if all(~isinf(phi.spaceRadius))
	qStartInt = ceil(phi.spaceRadius(1)/yUpStep);
	qSize = 2*qStartInt + 1;
	qStep = yUpStep;
	qStart = -qStartInt * yUpStep;
else
	qSize = 2*yUpSize - 1;
	qStep = yUpStep;
	qStart = -(yUpSize-1) .* yUpStep;
end
qGrid = makeGrid(qStart, qStep, qSize);


% interp setup
HTg = zeros(1, prod(xSize));
xkGrid = makeGrid(xStart, xStep, xSize); % might be big! :(
Binterp = BSplineKernel(degree, D-1, 1, 1);
gUp = zeros( [yUpSize 1] );

%% main loop
% (we chose to process projInd by projInd instead of doing the upsamping and
% convolution in a single step to reduce memory consumption)
for projInd = 1:numProj
	
	% project the xGrid
	Px = Ps(:, :, projInd) * xkGrid;
	
	% upsample g
	gUp( upSampleVects{:} ) = g( yInds{:}, projInd );
	
	% sample the x-ray transform of the kernel
	q = phi.xray(qGrid, Ps(:,:, projInd));
	q = reshape(q, [qSize 1]);
	
	% compute g * phi
	if D == 2
		gConvPhi = conv(gUp, q, 'full');
	elseif D == 3
		
		gConvPhi = conv2(gUp, q, 'full');
	end
	
	% interpolate	
	c = gConvPhi;
	if D==2
		c = coeff( c', degree);
	elseif D==3 % separable spline fitting
		c = spline_fit(c, degree);
		c = spline_fit( c', degree)';
	end % todo: for some reason, coeff is not sufficiently accurate for 2D interpolation and the fit _at the points of q_ is only 60 dB.
	
	if D == 2
		vals = Binterp.reconstruct(c, Px/yUpStep - (yUpGrid(1)/yUpStep -(length(q)-1)/2), 'mirror');
	elseif D == 3
		vals = Binterp.reconstruct(c, bsxfun(@minus, bsxfun(@times, Px, 1./yUpStep'),  yUpStart'./yUpStep' - (size(q)'-1)/2 ), 'mirror');
	end
	
	HTg = HTg + vals;
end
HTg = reshape(HTg, xSize);

global debugFlag
if debugFlag
	F = @(x) Binterp.reconstruct(c, x/yUpStep - (yUpGrid(1)/yUpStep -(length(q)-1)/2), 'mirror');
	x = sort([(Px(:)); linspace(min(Px(:)), max(Px(:)), length(Px))'; yUpGrid'])';
	h = plot(x, F(x));
	plot(yUpGrid, F(yUpGrid), 'o', 'color', h.Color);
	set(gca, 'colororderIndex', get(gca, 'colorOrderIndex')-1);
	
end
