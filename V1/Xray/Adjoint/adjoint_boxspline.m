function HTg = adjoint_boxspline(g, yStart, yStep, ySize, Ps, ...
    xStart, xStep, xSize, box)
xkGrid = makeGrid(xStart, xStep, xSize);
ymGrid = makeGrid(yStart, yStep, ySize);
[allshift,allcoe]= box.func(Ps);
if ndims(Ps) == 3
    ps(:,:)=-Ps(1,:,:);
else
    ps(:,:)=-Ps(1,:,:)';
end
offset = ps' * xkGrid;

HTg=boxsplineinterp_reduce(allshift,allcoe,g,offset',ymGrid,(xStep*abs(ps)+box.acqwidth)/2,yStep,box.acqwidth);
HTg=HTg./xStep(1)/xStep(1);

