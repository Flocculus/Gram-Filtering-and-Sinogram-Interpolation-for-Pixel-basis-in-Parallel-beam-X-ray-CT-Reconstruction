classdef EllipsoidPhantom < XRayPhantom
    %modified on 3/31/2019 add conv(projection,boxspline)
    properties
        E;
    end
    
    methods(Static)
        function f = example(num)
            % example of the usage of the ellipsoid phantom
            if ~exist('num', 'var') || isempty(num)
                num = 1;
            end
            
            switch num
                %----------------------------------------------------------------
                case 1
                    % each row is an ellipse.
                    % intensity, a, b, centerX, centerY, rotation, just like 'help
                    % phantom'
                    ells = [1, 10, 3, 11, 5, pi/8
                        .4, 5, 5, 0, 0, 0
                        2, 10, 7, -1, -10, 2*pi/3
                        2 .5 .5 -5 7 0]
                    p = EllipsoidPhantom(ells);
                    
                    xStart = [-25 -25];
                    xStep = [.25 .25];
                    xSize = [200 200];
                    f = p.eval(xStart, xStep, xSize);
                    
                    figure
                    % rows of f correspond to x_0, cols to x_1 , pages to x_2, so we need to
                    % transpose to plot it in the conventional way, with "x" going
                    % left-right and "y" going up-down
                    imagesc( (0:xSize(1)-1) * xStep(1) + xStart(1), ...
                        (0:xSize(2)-1) * xStep(2) + xStart(2), ...
                        f');
                    axis xy
                    colorbar
                    xlabel('x_0')
                    ylabel('x_1')
                    %----------------------------------------------------------------
                case 2 % x-ray projections of a simple 2D phantom
                    ells = [1, 10, 3, 11, 5, pi/10 + pi/2]
                    p = EllipsoidPhantom(ells);
                    
                    
                    %% calculate a sinogram
                    numProj = 10; % number of projection angles
                    P = zeros(1, 2, numProj);
                    t = (0:numProj-1) * pi/numProj;
                    for i = 1:numProj
                        P(:,:,i) = [cos(t(i)) sin(t(i))];
                    end
                    
                    yStart = -25;
                    yStep = 1;
                    ySize = 50;
                    g = p.xray(yStart, yStep, ySize, P);
                    
                    sino_h = figure;
                    imagesc(g, 'XData', t, 'YData', (0:ySize-1)*yStep + yStart)
                    xlabel('perpendicular to projection angle')
                    ylabel('y')
                    axis manual;
                    title('sinogram');
                    
                    %% display the phanotm
                    xStart = [-25 -25];
                    xStep = [.25 .25];
                    xSize = [200 200];
                    f = p.eval(xStart, xStep, xSize);
                    
                    f_h = figure;
                    % rows of f correspond to x_0, cols to x_1 , pages to x_2, so we need to
                    % transpose to plot it in the conventional way, with "x" going
                    % left-right and "y" going up-down
                    imagesc((0:xSize(1)-1) * xStep(1) + xStart(1), ...
                        (0:xSize(2)-1) * xStep(2) + xStart(2), ...
                        f');
                    axis xy
                    colorbar
                    axis manual % so the axes don't change when we plot lines
                    xlabel('x_0')
                    ylabel('x_1')
                    
                    %% display some lines
                    
                    styles = {'g', 'y', 'r'};
                    
                    inds = [2 7 7
                        38 18 37];
                    for i = 1:size(inds,2)
                        pInd = inds(1,i);
                        yInd = inds(2,i);
                        
                        figure(f_h)
                        hold on
                        y = makeGrid(yStart, yStep, ySize);
                        x = P(:,:,pInd)' * y;
                        u = null( P(:,:,pInd) );
                        %plot(x(1,yInd), x(2,yInd), 'o')
                        l = [x(:,yInd) + xStep(1)*xSize(1)*u, x(:,yInd) - xStep(2)*xSize(2)*u];
                        plot(l(1,:), l(2,:), styles{i});
                        
                        figure(sino_h)
                        hold on
                        plot(t(pInd), y(yInd), ['.' styles{i}], 'markersize', 20);
                    end
                    
                    
                    %----------------------------------------------------------------
                case 3 % 3-d example
                    % each row is an ellipse.
                    ells = [1, 5, 3, 9, 0, 0, 0, 0, 0, 0
                        1, 7, 17, 3, 10, 3, 0, 0, 0, 0
                        3, 10, 1, 1, -10, -5, 0, pi/4, 0, 0
                        2, 10, 4, 20, -15, 10, 0, pi/8, 2*pi/3, 0]
                    
                    p = EllipsoidPhantom(ells);
                    
                    xStart = [-25 -25 -25];
                    xStep = [.25 .25 .25];
                    xSize = [200 200 200];
                    f = p.eval(xStart, xStep, xSize);
                    
                    
                    for slice = [65, 101, 121];
                        figure
                        imagesc((0:xSize(1)-1) * xStep(1) + xStart(1), ...
                            (0:xSize(2)-1) * xStep(2) + xStart(2), ...
                            f(:,:,slice)', [0 5]);
                        axis xy
                        colorbar
                        xlabel('x_0')
                        ylabel('x_1')
                        title(sprintf('x_2 = %g', xStep(3)*(slice-1) + xStart(3)));
                    end
                    
                    %% take x-ray measurements
                    % setup the y grid, which is 2D
                    yStart = [-30 -30];
                    yStep = [.25 .25];
                    ySize = [300 300];
                    
                    % set angles determining the projection directions
                    phi =   [0  0     0     pi/4  pi/2  pi/4];
                    theta = [0  -pi/2  pi/2  0     0     pi/4];
                    psi =   [0  0     0     0     0     0   ];
                    numProj = length(phi); % number of projection angles
                    
                    % generate a basis orthogonal to the projection direction for
                    % each projection.
                    [r0, r1, r2] = Euler3D(phi, theta, psi);
                    P = reshape([r0(:)'; r1(:)'], 2, 3, numProj);
                    
                    % take the measurements
                    g = p.xray(yStart, yStep, ySize, P);
                    
                    % plot
                    
                    
                    for gInd = [1 2];
                        figure
                        imagesc( (0:ySize(1)-1)*yStep(1) + yStart(1), ...
                            (0:ySize(2)-1)*yStep(2) + yStart(2), ...
                            g(:,:,gInd)');
                        axis xy
                        colorbar
                        xlabel(sprintf('y_0 = [%.0f %.0f %.0f]', r0(:,gInd)));
                        ylabel(sprintf('y_1 = [%.0f %.0f %.0f]', r1(:,gInd)));
                        title(sprintf('projection along [%.0f %.0f %.0f]', r2(:,gInd)));
                    end
                    
                    figure
                    volume_plot(xStart, xStep, xSize, f, .5);
                    
                case 4
                    % 3D shepp phantom from "Computerized tomography and nuclear
                    % magnetic resonance" 1980
                    
                    % copied directly from the paper
                    sheppFull = [
                        0 0 0 .7233 .9644 1.27 2 1 0 0 0 1 0 0 0 1
                        0 -.0184 -.0185 .7008 .9246 1.2241 -.98 1 0 0 0 1 0 0 0 1
                        .2583 .7534 0 .127 .127 .127 -1 1 0 0 0 1 0 0 0 1
                        -.2583 .7534 0 .127 .127 .127 -1 1 0 0 0 1 0 0 0 1
                        0 1.1398 -.1957 .1270 .34 .17 1.5 1 0 0 0 .5446 -.8387 0 .8387 .5446
                        0 0 -.7620 .4575 .6099 .5080 -1 1 0 0 0 1 0 0 0 1
                        .7076 -.1378 -.1905 .0635 .3175 .3175 1 .9903 -.1085 -.0865 .1089 .9941 0 .086 -.0094 .9963
                        -.7076 -.1378 -.1905 .0635 .3175 .3175 1 -.9903 -.1085 -.0865 -.1089 .9941 0 -.086 -.0094 .9963
                        -.08 -.6050 .3810 .046 .023 .023 .01 1 0 0 0 1 0 0 0 1
                        0 -.6050 .3810 .023 .023 .046 .01 1 0 0 0 1 0 0 0 1
                        .06 -.605 .381 .023 .046 .046 .01 1 0 0 0 1 0 0 0 1
                        0 .1 .381 .046 .046 .046 .01 1 0 0 0 1 0 0 0 1
                        0 -.1 .127 .2581 .2581 .2581 .01 1 0 0 0 1 0 0 0 1
                        0 .35 .381 .21 .25 .23 .01 1 0 0 0 1 0 0 0 1
                        .22 0 .381 .110 .31 .254 -.02 .9511 -.3090 0 .3090 .9511 0 0 0 1
                        -.22 0 .381 .1600 .4100 .3810 -.02 -.9511 -.3090 0 -.3090 .9511 0 0 0 1
                        .5600 -.4000 .3810 .03 .2 .2 .03 .9192 -.3381 .2020 .3452 .9385 0 .1896 -.0697 -.9794
                        ];
                    shepp = [
                        %0 0 0 .7233 .9644 1.27 2 1 0 0 0 1 0 0 0 1
                        %0 -.0184 -.0185 .7008 .9246 1.2241 -.98 1 0 0 0 1 0 0 0 1
                        %.2583 .7534 0 .127 .127 .127 -1 1 0 0 0 1 0 0 0 1
                        %-.2583 .7534 0 .127 .127 .127 -1 1 0 0 0 1 0 0 0 1
                        %0 1.1398 -.1957 .1270 .34 .17 1.5 1 0 0 0 .5446 -.8387 0 .8387 .5446
                        %0 0 -.7620 .4575 .6099 .5080 -1 1 0 0 0 1 0 0 0 1
                        %.7076 -.1378 -.1905 .0635 .3175 .3175 1 .9903 -.1085 -.0865 .1089 .9941 0 .086 -.0094 .9963
                        %-.7076 -.1378 -.1905 .0635 .3175 .3175 1 -.9903 -.1085 -.0865 -.1089 .9941 0 -.086 -.0094 .9963
                        -.08 -.6050 .3810 .046  .023 .023  .01 1 0 0 0 1 0 0 0 1
                        0    -.6050 .3810 .023  .023 .046  .01 1 0 0 0 1 0 0 0 1
                        .06   -.605 .381  .023  .046 .046  .01 1 0 0 0 1 0 0 0 1
                        0      .1   .381  .046  .046 .046  .01 1 0 0 0 1 0 0 0 1
                        0      -.1  .127 .2581  .2581 .2581 .01  1 0 0 0 1 0 0 0 1
                        0      .35  .381 .21    .25   .23   .01  1 0 0 0 1 0 0 0 1
                        .22    0    .381 .110 .31 .254      -.02 .9511 -.3090 0 .3090 .9511 0 0 0 1
                        -.22   0    .381 .1600 .4100 .3810   -.02 -.9511 -.3090 0 -.3090 .9511 0 0 0 1
                        .5600  -.4000 .3810 .03 .2 .2       .03   .9192 -.3381 .2020 .3452 .9385 0 .1896 -.0697 -.9794
                        ];
                    
                    
                    % massage into our format
                    ells = zeros(size(shepp, 1), 16);
                    ells(:, 1) = shepp(:, 7);
                    ells(:, 2:4) = shepp(:, 1:3);
                    ells(:, 5:7) = shepp(:, 4:6);
                    ells(:, 8:16) = shepp(:, 8:16);
                    
                    p = EllipsoidPhantom(ells);
                    
                    xStart = [-.7 -.7 -.16];
                    xStep = [.01 .01 .01];
                    xSize = [141 131 101];
                    f = p.eval(xStart, xStep, xSize);
                    
                    figure
                    prettyPlot(xStart, xStep, xSize, f);
                    
                    % set angles determining the projection directions
                    % fixed axis is x_2
                    phi = linspace(0, pi, 1201);
                    theta = -pi/2 * ones(size(phi));
                    psi =   zeros(size(phi));
                    numProj = length(phi); % number of projection angles
                    
                    % generate a basis orthogonal to the projection direction for
                    % each projection.
                    [r0, r1, r2] = Euler3D(phi, theta, psi);
                    Ps = reshape([r1(:)'; r0(:)'], 2, 3, numProj); % swapping r1 r0 on purpose
                    
                    % take the measurements
                    yStart = [-.8 -.16];
                    yStep = [.0025 .0025]; % y isn't symmetrical (-.8 + 600*.0025 =/= .8), why? - well it doesn't need to be bigger
                    ySize = [601, 401];
                    g = p.xray(yStart, yStep, ySize, Ps);
                    
                    save('shepp3d', 'p', 'g', 'yStart', 'yStep', 'ySize', 'Ps', '-v7.3');
                    
                    % add noise w/ SNR of ... 0db
                    targetSNR = -10;
                    n = 10^(-targetSNR/20) * randn(size(g)) * sqrt(sum(g(:).^2) / numel(g)) ;
                    g = g + n;
                    
                    save('shepp3d_noisey', 'p', 'g', 'yStart', 'yStep', 'ySize', 'Ps', '-v7.3');
                    
            end % switch
            
            
        end % example func
    end % static methods
    
    methods %----------------------------------------------------------------
        
        function obj = EllipsoidPhantom(E)
            % obj = EllipsoidPhantom(E)
            %
            % E - defines the ellipsiods in the phantom.
            %     FOR A 2D PHANTOM: E is P X 6, where P is the number of ellipses
            %     E(p, 1) = intensity
            %     E(p, 2) = semi-axis 1 length
            %     E(p, 3) = semi-axis 2 length
            %     E(p, 4) = center coordinate 1
            %     E(p, 5) = center coordinate 2
            %     E(p, 6) = angle (radians) between axis 1 and the horizontal
            %
            %     FOR A 3D PHANTOM: E is P X 10
            %     E(p, 1) = intensity
            %     E(p, 2) = semi-axis 1 length
            %     E(p, 3) = semi-axis 2 length
            %     E(p, 4) = semi-axis 3 length
            %     E(p, 5) = center coordinate 1
            %     E(p, 6) = center coordinate 2
            %     E(p, 7) = center coordinate 3
            %     E(p, 8:10) = phi, theta, psi giving the z-y-z euler angles for
            %     the semi-axes
            %
            %     Or alternately:
            %     E(p, 1) = intensity
            %     E(p, 2) = center coordinate 1
            %     E(p, 3) = center coordinate 2
            %     E(p, 4) = center coordinate 3
            %     E(p, 5:7) = semi-axis lengths
            %     E(p, 8:10) = semi-axis 1
            %     E(p, 11:13) = semi-axis 2
            %     E(p, 14:16) = semi-axis 3
            
            obj.E = E;
            if size(E, 2) == 6
                obj.D = 2;
            elseif size(E, 2) == 10 || size(E, 2) == 16
                obj.D = 3;
            else
                error('unknown format for input E');
            end
        end
        
        function f = eval(obj, xStart, xStep, xSize)
            % f = eval(obj, xStart, xStep, xSize)
            
            if length(xStart) ~= obj.D || ...
                    length(xStep) ~= obj.D || ...
                    length(xSize) ~= obj.D
                error('dimension mismatch');
            end
            
            x = makeGrid(xStart, xStep, xSize);
            f = reshape(obj.evalPoints(x), [xSize, 1]);
            
        end
        
        function f = evalPoints(obj, x)
            % setup output grid
            
            f = zeros(size(x,2), 1);
            
            % do calcs for each ellipse
            for eInd = 1:size(obj.E,1);
                % find the rotation matrix, Phi'
                if obj.D == 2
                    R = [cos(obj.E(eInd, 6)) -sin(obj.E(eInd, 6))
                        sin(obj.E(eInd, 6)) cos(obj.E(eInd, 6))]';
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 10
                    phi = obj.E(eInd, 1+2*obj.D+1);
                    theta = obj.E(eInd, 1+2*obj.D+2);
                    psi = obj.E(eInd, 1+2*obj.D+3);
                    [r1, r2, r3] = Euler3D(phi, theta, psi);
                    R = [r1 r2 r3]';
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 16
                    R = reshape(obj.E(eInd, 8:16), 3, 3)';
                    axes = obj.E(eInd, 5:7)';
                    offset = obj.E(eInd, 2:4)';
                end
                
                
                % handle the center offset
                pos = R * bsxfun(@minus, x, offset);
                
                % handle the rotation
                vals = sum( bsxfun(@times, pos.^2, 1./axes.^2) );
                
                % add ellipse to the output
                f(vals <= 1) = f(vals <= 1) + obj.E(eInd,1);
            end
            
            
            
        end
        
        function g = xray(obj, yStart, yStep, ySize, P)
            % g = xray(obj, yStart, yStep, ySize, P)
            
            if length(yStart) ~= obj.D-1 || ...
                    length(yStep) ~= obj.D-1 || ...
                    length(ySize) ~= obj.D-1
                error('dimension mismatch');
            end
            
            
            y = makeGrid(yStart, yStep, ySize);
            g = zeros([ySize size(P, 3)]);
            
            % for each ellipse...
            for eInd = 1:size(obj.E,1);
                % get the center and axis lengths
                offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                axes = obj.E(eInd, (1:obj.D) + 1)';
                
                % find the rotation matrix for the ellipse
                if obj.D == 2
                    R = [cos(obj.E(eInd, 6)) -sin(obj.E(eInd, 6))
                        sin(obj.E(eInd, 6)) cos(obj.E(eInd, 6))]';
                    
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 10
                    phi = obj.E(eInd, 1+2*obj.D+1);
                    theta = obj.E(eInd, 1+2*obj.D+2);
                    psi = obj.E(eInd, 1+2*obj.D+3);
                    [r1, r2, r3] = Euler3D(phi, theta, psi);
                    R = [r1 r2 r3]';
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 16
                    R = reshape(obj.E(eInd, 8:16), 3, 3)';
                    axes = obj.E(eInd, 5:7)';
                    offset = obj.E(eInd, 2:4)';
                end
                
                % for each projection angle
                for pInd = 1:size(P,3)
                    % find the intersection of the ellipse with the line
                    % y(t) = x + t*u
                    
                    % we do this in the coordinate system of the ellipse, so that it
                    % is axis-aligned.
                    x = R * bsxfun(@minus, P(:,:,pInd)' * y, offset);
                    u = R * null( P(:,:,pInd) );
                    
                    % find real roots of the quadratic defined by the intersection
                    % Ax^2 + Bx + C = 0
                    A = sum(u.^2 ./ axes.^2);
                    B = 2 * sum( bsxfun(@times, u./axes.^2, x) );
                    C = sum( bsxfun(@times, x.^2, 1 ./ axes.^2) ) - 1;
                    
                    d = sqrt( B.^2 - 4 * A * C );
                    r1 = ( -B + d )/ (2*A);
                    r2 = ( -B - d )/ (2*A);
                    
                    h = imag(r1) == 0 & imag(r2) == 0;
                    
                    hInd = find(h)' + prod(ySize) * (pInd-1)';
                    g(hInd) = g(hInd) + abs(r1(h) - r2(h))' * obj.E(eInd,1);
                end
            end
            
        end
        
        function g = convoledxray(obj, yStart, yStep, ySize, P,step)
            % g = xray(obj, yStart, yStep, ySize, P)
            
            if length(yStart) ~= obj.D-1 || ...
                    length(yStep) ~= obj.D-1 || ...
                    length(ySize) ~= obj.D-1
                error('dimension mismatch');
            end
            
            
            y = makeGrid(yStart, yStep, ySize);
            g = zeros([ySize size(P, 3)]);
            
            % for each ellipse...
            for eInd = 1:size(obj.E,1)
                eInd
                %tic
                % get the center and axis lengths
                offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                axes = obj.E(eInd, (1:obj.D) + 1)';
                
                % find the rotation matrix for the ellipse
                if obj.D == 2
                    R = [cos(obj.E(eInd, 6)) -sin(obj.E(eInd, 6))
                        sin(obj.E(eInd, 6)) cos(obj.E(eInd, 6))]';
                    
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 10
                    phi = obj.E(eInd, 1+2*obj.D+1);
                    theta = obj.E(eInd, 1+2*obj.D+2);
                    psi = obj.E(eInd, 1+2*obj.D+3);
                    [r1, r2, r3] = Euler3D(phi, theta, psi);
                    R = [r1 r2 r3]';
                    
                    axes = obj.E(eInd, (1:obj.D) + 1)';
                    
                    offset = obj.E(eInd, (1:obj.D)+ obj.D+1)';
                    
                elseif size(obj.E, 2) == 16
                    R = reshape(obj.E(eInd, 8:16), 3, 3)';
                    axes = obj.E(eInd, 5:7)';
                    offset = obj.E(eInd, 2:4)';
                end
                
                % for each projection angle
                for pInd = 1:size(P,3)
                    pInd
                    % find the intersection of the ellipse with the line
                    % y(t) = x + t*u
                    
                    % we do this in the coordinate system of the ellipse, so that it
                    % is axis-aligned.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%numericial
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%method
%                                         x = @(y) (R * (P(:,:,pInd)' * y - offset));
%                                         u = R * null( P(:,:,pInd) );
%                                         % find real roots of the quadratic defined by the intersection
%                                         % Ax^2 + Bx + C = 0
%                                         A = sum(u.^2 ./ axes.^2);
%                                         B =@(y) 2 * sum( (u./axes.^2).*x(y) );
%                                         C =@(y) sum( (x(y).^2).* (1 ./ axes.^2))  - 1;
%                                        
%                                         f = @(x) (sqrt( B(x).^2 - 4 * A * C(x) ) / A);
%                                         start = y(1);
%                                         b = @(yy) yy;
%                                         a = ones(1,length(b))*(start);
%                                         g1 = @(x,yy) (b(yy+step/2)-a).*f(a+(b(yy+step/2)-a)*x);
%                                         g2 = @(x,yy) (b(yy-step/2)-a).*f(a+(b(yy-step/2)-a)*x);
%                                         
%                                         I_integralv =@(yy) integral(@(x) g1(x,yy)-g2(x,yy), 0, 1, 'ArrayValued', true);
%                                         sol = I_integralv(y);
                    %                     sol
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%symbolic
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%method
                    %
                    %
                    %                                         syms f(t) s
                    %                                         u = R * null( P(:,:,pInd) );
                    %                                          A = sum(u.^2 ./ axes.^2);
                    %                                         f(t) =  (sqrt( (2 * sum( (u./axes.^2).*(R * (P(:,:,pInd)' * t - offset)) )).^2 - 4 * A * (sum( ((R * (P(:,:,pInd)' * t - offset)).^2).* (1 ./ axes.^2))  - 1) ) / A);
                    %                                         start = y(1);
                    %                                         fff = int(f,[start s+0.5*step])-int(f,[start s-0.5*step]);
                    %                                         s=y;
                    %                                         sol = eval(fff);
                    %
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% put
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% offset in
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to the
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expression
%                     
%                     syms f(t) s cosA sinA
%                     
%                     u = (R * [-sinA; cosA] );
%                     A = sum(u.^2 ./ axes.^2);
%                     f(t) =  simplify((sqrt( (2 * sum( (u./axes.^2).*(R * ([cosA;sinA] * t - offset)) )).^2 - 4 * A * (sum( ((R * ([cosA;sinA] * t - offset)).^2).* (1 ./ axes.^2))  - 1) ) / A));
%                     cosA = P(1,1,pInd);
%                     sinA = P(1,2,pInd);
%                     start = y(1);
%                     fff = eval(int(simplify(eval(f(t))),[start s+0.5*step])-int(simplify(eval(f(t))),[start s-0.5*step]));
%                     s=y;
%                     sol = eval(fff);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% take out
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% the
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% offset
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    syms f(t) s cosA sinA
                    
                    u = (R * [-sinA; cosA] );
                    A = sum(u.^2 ./ axes.^2);
                    f(t) =  simplify((sqrt( (2 * sum( (u./axes.^2).*(R * ([cosA;sinA] * t )) )).^2 - 4 * A * (sum( ((R * ([cosA;sinA] * t )).^2).* (1 ./ axes.^2))  - 1) ) / A));
                    cosA = P(1,1,pInd);
                    sinA = P(1,2,pInd);
                    start = y(1);
                    fff = eval(int(simplify(eval(f(t))),[start s+0.5*step])-int(simplify(eval(f(t))),[start s-0.5*step]));
                    %fff = eval(int(simplify(eval(f(t))),[s-0.5*step s+0.5*step]));
                    shift = offset'*[cosA; sinA];
                    
                    s=y-shift;
                    sol = eval(fff);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    h = real(sol) ~= 0;
                    hInd = find(h)' + prod(ySize) * (pInd-1)';
                    g(hInd) = g(hInd) + real(sol(h))' * obj.E(eInd,1)/step;
                    
                    % h = imag(r1) == 0 & imag(r2) == 0;
                    
                    %hInd = find(h)' + prod(ySize) * (pInd-1)';
                    % g(hInd) = g(hInd) + abs(r1(h) - r2(h))' * obj.E(eInd,1);
                end
                %toc
            end
            
        end
    end
    
    
end