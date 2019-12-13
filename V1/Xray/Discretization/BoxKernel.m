classdef BoxKernel
    properties
        xStep;
        acqwidth;
    end
    
    
    methods %----------------------------------------------------------------
        
        function obj = BoxKernel(xStep,width)
            obj.xStep = xStep;
            obj.acqwidth = width;
        end
        
        function g = convoledxray(obj, args,P)
            Projection=BoxSpline1D3d(P,obj.xStep(1),obj.acqwidth);
            g = Projection{1}(args);
        end

        function [shift,coe] = func_reduce(obj,P)
            
            step = obj.xStep(1);
            width = obj.acqwidth;
            if (ndims(P) ~= 2)
                theta(:,:) = abs(P(1,:,:));
                xi1 = theta(1,:);
                xi2 = theta(2,:);
            else
                theta(:,:) = abs(P);
                xi1 = theta(1);
                xi2 = theta(2);
            end
            
            
            
            n=size(xi1,2);
            shift=cell(n,1);
            coe=cell(n,1);
            if width ~= 0
                for i=1:n
                    tmp = [xi1(i) xi2(i) width];
                    if ((xi1(i)) == 1 || (xi2(i)) == 1)%3d,width(implement in boxspline_reduce.c)|width|step
                        shifttmp = [0 step width step+width]';
                        [shift{i},I] = sort(shifttmp);
                        coetmp=(step/width)*[1 -1 -1 1];
                        coe{i} = coetmp(I);
                    else
                        para = abs(1/tmp(1)/tmp(2));%4d,width(implement in boxspline_reduce.c)|width|cos,sin
                        shifttmp = [0 tmp(1)*step tmp(2)*step width tmp(1)*step+tmp(2)*step tmp(1)*step+width tmp(2)*step+width tmp(1)*step+tmp(2)*step+width]';
                        [shift{i},I] = sort(shifttmp);
                        coetmp = (1/width)*para*[1 -1 -1 -1 1 1 1 -1];
                        coe{i} = coetmp(I);
                        %shift{i} = [0 tmp(1)*step tmp(2)*step width tmp(1)*step+tmp(2)*step tmp(1)*step+width tmp(2)*step+width tmp(1)*step+tmp(2)*step+width]'-(tmp(1)*step+tmp(2)*step+width)/2;
                        %coe{i} = (step/width)*para*[1 -1 -1 -1 1 1 1 -1]'/step;
                    end
                end
            else
                for i=1:n
                    tmp = [xi1(i) xi2(i)];%2d,width(implement in boxspline_reduce.c)|step
                    if ((xi1(i)) == 1 || (xi2(i)) == 1)
                        shift{i} = [0 step]';
                        coe{i} = [step -step];
                    else
                        para = abs(1/tmp(1)/tmp(2));%3d,width(implement in boxspline_reduce.c)|cos,sin
                        shifttmp = [0 tmp(1)*step tmp(2)*step tmp(1)*step+tmp(2)*step]';
                        [shift{i},I] = sort(shifttmp);
                        %coe(i) = para;
                        coetmp=(para)*[1 -1 -1 1];
                        coe{i} = coetmp(I);
                        %shift{i} = [0 tmp(1)*step tmp(2)*step width tmp(1)*step+tmp(2)*step tmp(1)*step+width tmp(2)*step+width tmp(1)*step+tmp(2)*step+width]'-(tmp(1)*step+tmp(2)*step+width)/2;
                        %coe{i} = (step/width)*para*[1 -1 -1 -1 1 1 1 -1]'/step;
                    end
                end
            end
        end
        
        function [shift,coe] = func_exact(obj,P,yStep)
            
            step = obj.xStep(1);
            width = obj.acqwidth;
            if (ndims(P) ~= 2)
                theta(:,:) = abs(P(1,:,:));
                xi1 = theta(1,:);
                xi2 = theta(2,:);
            else
                theta(:,:) = abs(P);
                xi1 = theta(1);
                xi2 = theta(2);
            end
            
            
            
            n=size(xi1,2);
            shift=cell(n,1);
            coe=cell(n,1);
            if width ~= 0
                for i=1:n
                    tmp = [xi1(i) xi2(i) width];
                    if ((xi1(i)) == 1 || (xi2(i)) == 1)%4d,width,width|width,step,
                        shifttmp = [0 width 2*width 3*width step width+step 2*width+step 3*width+step]' - (3*width+step)/2;
                        [shift{i},I] = sort(shifttmp);
                        coetmp=(step/width/width)*[1 -3 3 -1 -1 3 -3 1];
                        coe{i} = coetmp(I);
                    else
                        para = abs(1/tmp(1)/tmp(2));%6d,width,width,width|width,cos,sin,
                        w = width;
                        c = tmp(1)*step;
                        s = tmp(2)*step;
                        shifttmp = [0 w 2*w 3*w 4*w c w+c 2*w+c 3*w+c 4*w+c s w+s 2*w+s 3*w+s 4*w+s s+c w+c+s 2*w+c+s 3*w+c+s 4*w+c+s ]' - (4*w+c+s)/2;
                        [shift{i},I] = sort(shifttmp);
                        coetmp = (1/width/width/width)*para*[1 -4 6 -4 1 -1 4 -6 4 -1 -1 4 -6 4 -1 1 -4 6 -4 1];
                        coe{i} = coetmp(I);
                    end
                end
            else
                for i=1:n
                    tmp = [xi1(i) xi2(i)];
                    if ((xi1(i)) == 1 || (xi2(i)) == 1)%2d,yStep|step
                        shifttmp = [0 yStep step yStep+step]' - (yStep+step)/2;
                        [shift{i},I] = sort(shifttmp);
                        %coetmp = (step/yStep)*[1 -1 -1 1];
                        coetmp = (step)*[1 -1 -1 1];
                        coe{i} = coetmp(I);
                    else
                        para = abs(1/tmp(1)/tmp(2));%4d,yStep,yStep|cos,sin
                        w = yStep;
                        c = tmp(1)*step;
                        s = tmp(2)*step;
                        shifttmp = [0 w 2*w c w+c 2*w+c s w+s 2*w+s c+s w+c+s 2*w+c+s]' - (2*w+c+s)/2;
                        [shift{i},I] = sort(shifttmp);
                        %coetmp=(para)/yStep/yStep*[1 -2 1 -1 2 -1 -1 2 -1 1 -2 1];
                        coetmp=(para)/yStep*[1 -2 1 -1 2 -1 -1 2 -1 1 -2 1];
                        coe{i} = coetmp(I);
                    end
                end
            end
        end
    end
    
    
end