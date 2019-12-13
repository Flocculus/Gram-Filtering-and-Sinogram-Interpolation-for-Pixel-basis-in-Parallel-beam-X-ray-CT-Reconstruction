function HTg = adjoint_boxspline_exact(g, yStart, yStep, ySize, Ps, ...
    xStart, xStep, xSize, box)
xkGrid = makeGrid(xStart, xStep, xSize);

[allshift,allcoe]= box.func_exact(Ps,yStep);
if ndims(Ps) == 3
    ps(:,:)=-Ps(1,:,:);
else
    ps(:,:)=-Ps(1,:,:)';
end
offset = ps' * xkGrid;

kernelwidth = xStep*abs(ps) + yStep + yStep + box.acqwidth + box.acqwidth;
kernelwidth(abs(kernelwidth - xStep(1) - yStep - yStep - box.acqwidth - box.acqwidth)<=1e-14) = kernelwidth(abs(kernelwidth - xStep(1) - yStep - yStep - box.acqwidth - box.acqwidth)<=1e-14) - yStep;
kernelwidth = kernelwidth/2;

if box.acqwidth ~= 0
    % width=0:
    % xStep*abs(ps)+ystep,
    % xStep*abs(ps)+ystep+ystep,
    % width~=0:
    % xStep*abs(ps)+acqwidth+acqwidth+ystep,
    % xStep*abs(ps)+acqwidth+acqwidth+ystep+ystep,
    [num_data,num_proj]=size(g);
    %Filter Q
    Q = sqrt(2)*(2*sqrt(2)-3).^(abs([-40:40]));
    
    %convolution
    filtered_g = zeros(num_data+80,num_proj);
    filtered_g(41:end-40,:)=g;
    %testg=filtered_g;
    for i = 1:num_proj
        if (abs(ps(1,i)) ~= 1 && abs(ps(2,i)) ~= 1)
            %filtered_g(:,i)=convxh(filtered_g(:,i),Q');
            filtered_g(:,i) = conv(g(:,i),Q);
%             testg(:,i)=convxh(testg(:,i),Q');
%             F1=@() conv(filtered_g(:,i),Q,'same');
%             timeit(F1)
%             F2=@() convxh(testg(:,i),Q');
%             timeit(F2)
        end
    end
    g = filtered_g;
    ySize = ySize+80;
    yStart = yStart - 40 * yStep;
end




ymGrid = makeGrid(yStart, yStep, ySize);


HTg=boxsplineinterp_exact(allshift,allcoe,g,offset',ymGrid,kernelwidth,yStep,box.acqwidth);
HTg=HTg./xStep(1)/xStep(1);

