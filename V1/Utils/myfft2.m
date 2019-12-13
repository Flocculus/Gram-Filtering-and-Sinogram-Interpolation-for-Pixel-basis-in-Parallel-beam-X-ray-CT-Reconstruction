function output = myfft2(A, b)
fftSize = size(A) + size(b) - 1;
h = zeros([fftSize 1]);
inds = make_grid_vectors(ones(1, length(size(A))) - floor(size(A)/2)-1, ones(1, length(size(A)), 1), size(A));
for i = 1:length(size(A))
    inds{i} = mod(inds{i}, size(h,i)) + 1;
end
h(inds{:}) = A;
output = fft2(h, fftSize(1), fftSize(2));