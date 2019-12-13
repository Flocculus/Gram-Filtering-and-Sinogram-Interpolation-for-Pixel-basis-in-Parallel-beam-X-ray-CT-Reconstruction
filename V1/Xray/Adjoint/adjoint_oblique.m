function HTg = adjoint_oblique(g, yStart, yStep, ySize, Ps, ...
	xStart, xStep, xSize, phi, yUpsampleRate, degree, oblique)
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


upSampleVects = make_grid_vectors(ones(1,D-1), yUpsampleRate, ySize);
yInds = make_grid_vectors(ones(1,D-1), ones(1,D-1), ySize);

interpDegree = degree;
if oblique
	preDegree = 0;
else
	preDegree = degree;
end

% find q, which is inner products of phi with the chosen B-spline
Bpre = BSplineKernel(preDegree, 1, 1, 1);
if all(~isinf(phi.spaceRadius))
	qStartInt = ceil((Bpre.spaceRadius*yUpStep + phi.spaceRadius(1))/yUpStep);
	qSize = 2*qStartInt + 1;
	qStep = yUpStep;
	qStart = -qStartInt * yUpStep;
else
	qSize = 2*yUpSize - 1;
	qStep = yUpStep;
	qStart = -(yUpSize-1) .* yUpStep;
end

qGridVects = make_grid_vectors(qStart, qStep, qSize);
points1 = qGridVects{1}(1:ceil(end/2)+1) - yUpStep(1)/2;

if D == 3
	points2 = qGridVects{2}(1:ceil(end/2)+1) - yUpStep(2)/2;
end

N=8;
kSize = round(qSize*sqrt(2)*N);
kStep = qStep/N;
kStart = -kSize.*kStep/2;
phiInterpDegree = 5;

if D == 2
	points = kStart + (0:kSize-1)*kStep;
	phiInterp = BSplineKernel(phiInterpDegree, 1, 1, qStep);
	xrayPhi = phi.xrayTransformFunc([]);
	q = xrayPhi.integrate(points); %this for 1D
	qc = coeff(q, phiInterpDegree);
elseif D == 3
	qc = cell(2,1);
	phiInterps = cell(2,1);
	for d = 1:D-1
		points = kStart(d) + (0:kSize(d)-1)*kStep(d);
		phiInterps{d} = BSplineKernel(phiInterpDegree, 1, 1, qStep(d));
		xrayPhi = phi.xrayTransformFunc([]);
		q = xrayPhi.integrate1D(points, d);
		qc{d} = coeff(q, phiInterpDegree);
	end
end



m0 = xrayPhi.m;
c0 = xrayPhi.c;


% interp setup
%HTg = zeros(xSize);
HTg = zeros(1, prod(xSize));
xkGrid = makeGrid(xStart, xStep, xSize); % might be big! :(


%% main loop
Binterp = BSplineKernel(interpDegree, 1, 1, 1);

gUp = zeros( [yUpSize 1] );
% (we chose to process projInd by projInd instead of doing the upsamping and
% convolution in a single step to reduce memory consumption)
for projInd = 1:numProj
	
	% project the xGrid
	Px = Ps(:, :, projInd) * xkGrid;
	
	
	if (isa(phi, 'IsotropicXRayKernel') && projInd==1) || ~isa(phi, 'IsotropicXRayKernel')
		xrayPhi = phi.xrayTransformFunc(Ps(:, :, projInd));
	end
	
	if preDegree == 0
		if D == 3
			ints = xrayPhi.integrate1D(points1',1) * xrayPhi.integrate1D(points2, 2);
			diffs = diff( diff(ints, 1, 1), 1, 2);
			q = 1/prod(yUpStep) * [diffs diffs(:, end-1:-1:1); diffs(end-1:-1:1, :) diffs(end-1:-1:1, end-1:-1:1)];		
			
			a = phiInterps{1}.reconstruct(qc{1}, (points1/(xrayPhi.m(1)/m0(1)) - kStart(1))/kStep(1), 'mirror');
			b = phiInterps{2}.reconstruct(qc{2}, (points2/(xrayPhi.m(2)/m0(2)) - kStart(2))/kStep(2), 'mirror');
			ints = a' * b;
			diffs = diff( diff(ints, 1, 1), 1, 2);
			q = 1/prod(yUpStep) * [diffs diffs(:, end-1:-1:1); diffs(end-1:-1:1, :) diffs(end-1:-1:1, end-1:-1:1)];		
		elseif D == 2
% 			ints = xrayPhi.integrate(points1);
% 			diffs =  diff(ints);
% 			q = 1/yUpStep * [diffs diffs(end-1:-1:1)]; % rely on symmetry
			
			%xrayPhi.c/c0
			ints =  phiInterp.reconstruct(qc, (points1/(xrayPhi.m/m0) - kStart)/kStep, 'mirror');
			diffs =  diff(ints);
			q = 1/yUpStep * [diffs diffs(end-1:-1:1)]; % rely on symmetry
			
		end
		
	else
		q = zeros(size(qGridVects{1}));
		for qInd = 1:length(qGridVects{1})
			q(qInd) = integral( @(x) Bpre.eval( (x - qGridVects{1}(qInd))/yUpStep ) .* xrayPhi.eval(x), ...
				qGridVects{1}(qInd)-yUpStep(1)*Bpre.spaceRadius, ...
				qGridVects{1}(qInd)+yUpStep(1)*Bpre.spaceRadius, ...
				'AbsTol', 1e-14, 'RelTol', 1e-10); % necessary to break 150 db
		end
		q = 1/yUpStep * q;
	end

	c = q;
	if D==2
		c = coeff( c', preDegree+interpDegree+1);
	elseif D==3 % separable spline fitting
		c = spline_fit(c, preDegree+interpDegree+1);
		c = spline_fit( c', preDegree+interpDegree+1)';
	end	% todo: for some reason, coeff is not sufficiently accurate for 2D interpolation and the fit _at the points of q_ is only 60 dB.
	
	% upsample g
	gUp( upSampleVects{:} ) = g( yInds{:}, projInd );
  
	% compute g conv phi
	if D == 2
		gConvc = conv(gUp, c, 'full');
	elseif D == 3
		gConvc = conv2(gUp, c, 'full');
	end
	% how to deal with y=0 pixel, where we should match exactly?
% 	tol = 1e-9;
% 	yUpGrid = makeGrid(yUpStart, yUpStep, yUpSize);
% 	isZero = all( abs(yUpGrid) < tol, 1);
% 	gConvc(isZero) = xrayPhi.eval(yUpGrid) * gUp(:);
	

	% interpolate
	
	if D == 2
		vals = Binterp.reconstruct(gConvc, Px/yUpStep - (yUpStart/yUpStep -(qSize-1)/2), 'mirror');
	elseif D == 3
		offset = yUpStart' ./ yUpStep' - (qSize'-1)/2;
		vals = Binterp.reconstruct(gConvc, bsxfun(@minus, bsxfun(@times, Px, 1./yUpStep'), ...
			offset ), 'mirror');
	end
	
	
	HTg = HTg + vals;
end
HTg = reshape(HTg, xSize);

global debugFlag
if debugFlag
	yUpGridVects = make_grid_vectors(yUpStart, yUpStep, yUpSize);
	%	plot(gConvPhiVects{1}, gConvPhi)
	F = @(x) Binterp.reconstruct(gConvc, x/yUpStep - (yUpGridVects{1}(1)/yUpStep -(qSize-1)/2), 'mirror');
	x = sort([(Px(:)); linspace(min(Px(:)), max(Px(:)), length(Px))'; yUpGridVects{1}'])';
	h = plot(x, F(x));
	plot(yUpGridVects{1}, F(yUpGridVects{1}), 'o', 'color', h.Color);
	set(gca, 'colororderIndex', get(gca, 'colorOrderIndex')-1);
	
end
