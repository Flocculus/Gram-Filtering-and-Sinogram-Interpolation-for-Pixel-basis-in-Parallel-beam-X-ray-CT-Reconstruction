function p = make_phantom(phantomInd, D)
rng(1)
make_enums;
offset = (1:D)'/3;
pDiameter = 8;
switch phantomInd
    case PHANTOMS.ellipse
        if D == 2
            %ells = [1, 2*( 1 + sqrt(5) ) / 4, 2.5, 1, 2, pi/4];
            ells = [1, 1, 1, 0, 0, 0];
        elseif D == 3
            ells = [1, 2*( 1 + sqrt(5) ) / 2, 2, 2,...
                0, 0, 0, pi/10 + pi/2, 0, 0];
        end
        p = EllipsoidPhantom(ells);
        
    case PHANTOMS.tinyElls
        ellLength = 1;
        numElls = 30;
        if D == 2
            ells = repmat([1, ellLength/8, ellLength/2, 0, 0, 0], numElls, 1);
            r = rand(numElls,1) * (pDiameter/2 - ellLength/2);
            t = rand(numElls,1)*2*pi;
            ells(:,4:5) = [r.*cos(t) r.*sin(t)];
            ells(:, 6) = rand(numElls, 1) * 2 * pi;
        elseif D == 3
            ells = repmat([1, ellLength/8, ellLength/8, ellLength/2, ...
                0, 0, 0, ... % center
                0, 0, 0], numElls, 1); % rotation
            r = rand(numElls,1) * (pDiameter/2 - ellLength/2);
            t = rand(numElls,2)*2*pi;
            ells(:,5:7) = [r.*cos(t) r.*sin(t)];
            ells(:, 6) = rand(numElls, 1) * 2 * pi;
            
        end
        p = EllipsoidPhantom(ells);
        
    case PHANTOMS.bigSinc
        p = GenFuncPhantom(IsotropicXRayKernel(SincFunc(), D, 2), offset);
        
    case PHANTOMS.box
        func = IsotropicGeneratingFunc(BSplineKernel(0), D, 1, 4);
        p = GenFuncPhantom(func, offset);
    case PHANTOMS.bigKBW
        func = IsotropicXRayKernel(KBWFunc(2,10.85,2), D, 1, 2);
        p = GenFuncPhantom(func, offset);
    case PHANTOMS.multiKBW
        spotAlpha = 10.85;
        spotRadius = 2;
        numSpots = 45;
        numColdSpots = 10;
        maxR = 1;
        minR = .35;
        func = cell(1,numSpots);
        for i = 1:numSpots
            func{i} = IsotropicXRayKernel(KBWFunc(2,spotAlpha,spotRadius), D, rand, rand*(maxR-minR) + minR);
        end
        for i = 1:numColdSpots
            func{i} = IsotropicXRayKernel(KBWFunc(2,spotAlpha,spotRadius), D, (2*rand)^2, minR);
        end
        
        rng(1);
        t = rand(1,numSpots)*2*pi;
        r = sqrt(rand(1,numSpots) * (pDiameter/2-maxR)^2);
        offset = [cos(t) .* r; sin(t) .* r] ;
        if D == 3
            offset(3,:) = 0;
        end
        p = GenFuncPhantom(func, offset);
    case PHANTOMS.smallKBW
        func1 = IsotropicXRayKernel(KBWFunc(2,10.85, 1), D, 1);
        func2 = IsotropicXRayKernel(KBWFunc(2,10.85, 1), D, -1);
        func3 = IsotropicXRayKernel(KBWFunc(2,10.85, 1), D, 1);
        p = GenFuncPhantom({func1, func2, func3}, [offset offset+1 offset+2]);
        
    case PHANTOMS.noise
        p = NoisePhantom(2);
    case PHANTOMS.smallSinc
        numSpots = 1;
        func = cell(1,numSpots);
        for i = 1:numSpots
            func{i} = SeparableSincKernel(D, ones(1,D)*.5,2*rand-1);
        end
        rng(1);
        offset = 2*(rand(D,numSpots)-.5) * 2;
        
        window = IsotropicXRayKernel(KBWFunc(2, 2.5, 4), D, 1);
        %p = GenFuncPhantom(window, [0;0]);
        %p = GenFuncPhantom(func, offset, window);
        p = GenFuncPhantom(func, [0;0]);
    case PHANTOMS.sheppSmooth
        load('sheppSmooth', 'vals', 'X', 'step');
        %X = X(:,65);
        %vals = vals(65);
        X = X*4;
        step = step * 4;
        if D == 3
            X(3,:) = 0;
        end
        
        funcs = cell(1,size(X,2));
        for i = 1:size(X,2);
            
            funcs{i} = IsotropicXRayKernel(KBWFunc(2, 10.85, 2), D, vals(i), step);
            %funcs{i} = IsotropicXRayKernel(KBWFunc(2, 10.85, 2), D, vals(i), 1);
            
        end
        
        p = GenFuncPhantom(funcs, X);
        
    case PHANTOMS.filament
        numSpots = 200;
        numSteps = 80;
        pathLength = 8;
        bias1 = ones(1,D)*.75;
        bias2 = ones(1,D)*.75;
        bias2(2) = -bias2(2);
        
        func = cell(1, numSpots*2);
        for i = 1:numSpots*2
            func{i} = IsotropicXRayKernel(KBWFunc(1,40, .6), D, 1);
        end
        
        steps = (bsxfun(@plus, randn(numSteps, D), bias1))/numSteps*pathLength;
        C = cumsum(steps);
        C = bsxfun(@minus, C, mean(C));
        
        Cinterp(:,1) = interp1((0:numSteps-1)/numSteps, C(:,1), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        Cinterp(:,2) = interp1((0:numSteps-1)/numSteps, C(:,2), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        if D == 3
            Cinterp(:,3) = interp1((0:numSteps-1)/numSteps, C(:,3), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        end
        
        steps = (bsxfun(@plus, randn(numSteps, D), bias2))/numSteps*pathLength;
        C = cumsum(steps);
        C = bsxfun(@minus, C, mean(C));
        
        Dinterp(:,1) = interp1((0:numSteps-1)/numSteps, C(:,1), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        Dinterp(:,2) = interp1((0:numSteps-1)/numSteps, C(:,2), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        if D == 3
            Dinterp(:,3) = interp1((0:numSteps-1)/numSteps, C(:,3), (0:numSpots-1)/numSpots, 'pchip', 'extrap');
        end
        
        
        p = GenFuncPhantom(func, [Cinterp' Dinterp']);
    case PHANTOMS.SheppLogan
        p = GenSheppLoganPhantom();
    case PHANTOMS.forbuild
        p = GenforbuildPhantom();
    case PHANTOMS.realImg
        p = ReadrealImg();
end
