t1 = clock;

outDir = 'Figures/';
addpath(genpath('Utils/'));
addpath(genpath('Invpblib/'));
addpath(genpath('FORBILDPhantom/'));
addpath(genpath('realimg/'));
addpath(genpath('Xray/'));
mex ./xray.c
mex ./Xray/Discretization/coeff.c
mex ./Xray/Discretization/interpol.c
mex ./Xray/Utils/ndgridMatrix.c
mex ./Xray/Adjoint/boxsplineinterp_reduce.c
mex ./Xray/Adjoint/boxsplineinterp_exact.c
mex ./grammatrix.c
mex ./convoledxray.c
mex ./realbackproject.c


rng(3);
D = 2;
numProjs = [15,30]; %number of projection directions
yStepMults = [1,2]; %Sampling step on sinogram domain compare to image domain. 0.5 means yStep = 0.5 xStep，it is also the width of the detector
xSizes =[ones(1,D) * 2;ones(1,D) * 4];%The resolution of the image
methods = [1 2 3 4 5]; % The methods{1:exact sinc 2:standard 3:orthogonal 4:oblique 5:box spline (with and without detector blur effect)}
phantomInds =[1,2]; %%{'ellipse', 'tinyElls', 'bigSinc', 'box', 'bigKBW','smallKBW', 'multiKBW', 'noise', 'smallSinc', 'sheppSmooth', 'filament', 'SheppLogan', 'forbuild', "realImg"}
resultbackprojection = zeros(length(numProjs),length(phantomInds),size(xSizes,1),2,length(yStepMults),length(methods),2);
%resultbackprojectionpic = zeros(length(numProjs),length(phantomInds),size(xSizes,1),2,length(yStepMults),length(methods),64,64);

%{
numProjects,
phantomIndex,
imgsize,
detector,[without, with]
ystepmults,
methods,
(backprojectionsnr,ssim),
%}
maxIter = 200;



for numProject = 1:length(numProjs)
    numProject
    thetaspace = [1:numProjs(numProject)]./numProjs(numProject)*180;
    %thetaspace = [60 30 0]
    thetas = thetaspace/180*pi;
    for phantomInd = 1:length(phantomInds)
        phantomInd
        for yStepMult = 1:length(yStepMults)
            yStepMult
            for xSize = 1:size(xSizes,1)
                 xSize
                
                %phantom
                pDiameter = 8; % phantom is nonzero inside the area
                p = make_phantom(phantomInds(phantomInd), D);
                % x grid
                xStep  = pDiameter ./ (xSizes(xSize,:));
                xStart = -pDiameter / 2 * ones(1,D) + xStep/2;
                % y grid
                yFudge = 2;
                yStep = xStep(1:D-1) * yStepMults(yStepMult);
                yStart = -pDiameter*yFudge/2+1e-10;
                ySize =  floor(pDiameter*yFudge./yStep) + 1;
                % projections
                Ps = zeros(1, 2, numProjs(numProject));
                for i = 1:numProjs(numProject)
                    Ps(:,:,i) = [cos(thetas(i)) sin(thetas(i))];
                end
                %get the opt discretization resultP12
                gtUpsampleRate = 6;
                xUpStep = xStep / gtUpsampleRate;
                xUpStart = -pDiameter / 2 * ones(1,D) + xUpStep/2;
                xUpSize = xSizes(xSize,:) * gtUpsampleRate;
                maxDiscIter = 50;
                c0 = zeros(xSizes(xSize,:));
                %boxspline
                H_Box = Sampler_Box(xStart, xStep, xSizes(xSize,:), gtUpsampleRate, pDiameter/2);

                gt = p.eval(xUpStart, xUpStep, xUpSize);
                %sinc
                phi = SeparableSincKernel(D, xStep);
                H_Sinc = Sampler(xStart, xStep, xSizes(xSize,:), gtUpsampleRate, phi, pDiameter/2);
                %get the forward projection data
                fprintf('start\n')
                tic
                g = xray(thetaspace,makeGrid(yStart, yStep, ySize),xUpSize(1),xUpStep(1),gt,yStep);
                gblur = convoledxray(thetaspace,makeGrid(yStart, yStep, ySize),xUpSize(1),xUpStep(1),gt,yStep);
                realbackprojection = realbackproject(thetaspace,xUpSize(1),xUpStep(1),gt,0,yStep);
                realbackprojectionblur = realbackproject(thetaspace,xUpSize(1),xUpStep(1),gt,1,yStep);
                toc
           
                    for method_ind = 1:size(methods,2)
                        method = methods(method_ind);
                        switch method
                            case 1 % exact
                                %get HTH
                                
                                G = XRay(xStart, xStep, xSizes(xSize,:), yStart, yStep, ySize, -Ps, phi, 'adjointMethod', method, 'adjointDegree', NaN, 'adjointUpsampleRate', NaN);
                                GTg = G'*g;
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,1) = snr(realbackprojection, realbackprojection - H_Sinc*GTg);
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,2)=ssim(H_Sinc*GTg,realbackprojection);
                                %resultbackprojectionpic(numProject,phantomInd,xSize,:,yStepMult,method,:,:) = GTg;
                            case 2 % standard interpolation
                                %get HTH
                          
                                G = XRay(xStart, xStep, xSizes(xSize,:), yStart, yStep, ySize, -Ps, phi, 'adjointMethod', method, 'adjointDegree', 3, 'adjointUpsampleRate', 1);
                                GTg = G'*g;
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,1) = snr(realbackprojection, realbackprojection - H_Sinc*GTg);
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,2)=ssim(H_Sinc*GTg,realbackprojection);
                                %resultbackprojectionpic(numProject,phantomInd,xSize,:,yStepMult,method,:,:) = GTg;
                            case 3 % orthogonal interpolation
                                %get HTH

                                G = XRay(xStart, xStep, xSizes(xSize,:), yStart, yStep, ySize, -Ps, phi, 'adjointMethod', method, 'adjointDegree', 3, 'adjointUpsampleRate', 1);
                                GTg = G'*g;
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,1) = snr(realbackprojection, realbackprojection - H_Sinc*GTg);
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,2)=ssim(H_Sinc*GTg,realbackprojection);
                                %resultbackprojectionpic(numProject,phantomInd,xSize,:,yStepMult,method,:,:) = GTg;
                            case 4 % oblique interpolation
                                %get HTH
                            
                                G = XRay(xStart, xStep, xSizes(xSize,:), yStart, yStep, ySize, -Ps, phi, 'adjointMethod', method, 'adjointDegree', 3, 'adjointUpsampleRate', 1);
                                GTg = G'*g;
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,1) = snr(realbackprojection, realbackprojection - H_Sinc*GTg);
                                resultbackprojection(numProject,phantomInd,xSize,:,yStepMult,method,2)=ssim(H_Sinc*GTg,realbackprojection);
                                %resultbackprojectionpic(numProject,phantomInd,xSize,:,yStepMult,method,:,:) = GTg;
                           case 5 %boxspine
                                for detectorInd = 1:2
                                    if (detectorInd == 1)%without detector blur
                                        g_real = g;
                                        box = BoxKernel(xStep,0);
                                        F = @() adjoint_boxspline_exact(g_real, yStart, yStep, ySize, Ps, xStart, xStep, xSizes(xSize,:), box);
                                        htg = F();
                                        resultbackprojection(numProject,phantomInd,xSize,detectorInd,yStepMult,method,1) = snr(realbackprojection, realbackprojection - H_Box*htg);
                                        resultbackprojection(numProject,phantomInd,xSize,detectorInd,yStepMult,method,2)=ssim(H_Box*htg,realbackprojection);
                                        %resultbackprojectionpic(numProject,phantomInd,xSize,detectorInd,yStepMult,method,:,:) = htg;
                                    else%with detector blur
                                        %get the delector blured forward projection data
                                        g_real = gblur;
                                        box = BoxKernel(xStep,yStep);
                                        F = @() adjoint_boxspline_exact(g_real, yStart, yStep, ySize, Ps, xStart, xStep, xSizes(xSize,:), box);
                                        htg = F();
                                        resultbackprojection(numProject,phantomInd,xSize,detectorInd,yStepMult,method,1) = snr(realbackprojectionblur, realbackprojectionblur - H_Box*htg);
                                        resultbackprojection(numProject,phantomInd,xSize,detectorInd,yStepMult,method,2)=ssim(H_Box*htg,realbackprojectionblur);
                                        %resultbackprojectionpic(numProject,phantomInd,xSize,detectorInd,yStepMult,method,:,:) = htg;
                                    end
                          
                                end %width
                        end%switch
                    end%method
            end %xSize
        end %yStepMults
    end %phantomInds
end %number of project
