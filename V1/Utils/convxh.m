function y = convxh(xIn, hIn, fftSize, isDFT)
% function y = convxh(xIn, hIn, fftSize, isDFT)
% compute y = h * x using zero boundary conditions and the fastest method I
% know (fft(ftt(...
%
% equivalent to y = convn(x, h, 'same'); but faster
%
% set isDFT true if you pass fft(h) as hIn, in which case fftSize =
% size(hIn)
% 
% the size is determined by the size of x
% the center of h is floor(size(h)/2) + 1
% works for x and h up to 3 dimensions
%
% also possible to pass hHat instead of h

% TODO: some callers may not deal well with this change in interface, but
% fixing them should be trivial.

if nargin < 3 || isempty(fftSize)
	fftSize = size(xIn) + size(hIn) - 1;
end
if nargin < 4 || isempty(isDFT)
	isDFT = false;
end

if isDFT
	hHat = hIn;
	fftSize = size(hHat);
else
	% handle the padding of hIn
	h = zeros([fftSize 1]);
	inds = make_grid_vectors(ones(1, length(size(hIn))) - floor(size(hIn)/2)-1, ones(1, length(size(hIn)), 1), size(hIn));
	for i = 1:length(size(hIn))
		inds{i} = mod(inds{i}, size(h,i)) + 1;
	end
	h(inds{:}) = hIn;
	
	if isvector(h)
		hHat  =  fft(h, max(fftSize));
	elseif ismatrix(h)
		hHat = fft2(h, fftSize(1), fftSize(2));
	elseif ndims(h) == 3
		hHat  =  fft(fft(fft(h,fftSize(1),1), fftSize(2),2), fftSize(3),3);
	else
		hHat =  fftn(h, fftSize);
	end

end


% do the convolution
if isvector(xIn)
	y  =  ifft(fft(xIn, max(fftSize)) .* hHat, 'symmetric');
elseif ismatrix(xIn)
	y = ifft2( fft2(xIn, fftSize(1), fftSize(2)) .* hHat, 'symmetric');
elseif ndims(xIn) == 3
	y  =  ifft(ifft(ifft( ...
		fft(fft(fft(xIn,fftSize(1),1), fftSize(2),2), fftSize(3),3) ...
		.* hHat ...
		,[],1), [],2), [],3,'symmetric');
else
	y = ifftn( fftn(xIn, fftSize) .* hHat, 'symmetric');
end

y = y(1:size(xIn,1), 1:size(xIn,2), 1:size(xIn,3));
%y = y(1:size(xIn,1), 1:size(xIn,2), 1:size(xIn,3));