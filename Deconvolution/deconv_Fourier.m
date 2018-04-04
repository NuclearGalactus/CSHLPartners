function output = deconv_Fourier(data, kernel, eps)
if nargin < 3
    eps = 0.1;
end
% data must be a row vector or a series of row vectors to be averaged
% kernel must be a single row vector
% eps must be a scalar, or absent to be defaulted to .1
% returns a row vector normalized so that its max value is 1
n = size(data,1); % length of data vector
L = size(data,2); % number of data vectors to process
Lx=size(data,2)-length(kernel)+1; % deconvolved length
Lx2=pow2(nextpow2(Lx)); % closed power of two to deconvolved length, makes fft faster
W = blackman(size(data,2)); % window with blackman filter to reduce spectral leakage
H= fft(kernel, Lx2); % transform kernel, kernel is constant so this only need be done once
out = zeros(n,Lx); % instantiate process vector
for counter = 1:n
 Y=fft(data(counter,:) .* W', Lx2); % transform data with window
 X = (Y.*conj(H))./(H.*conj(H)+eps*mean(H.*conj(H))); % deconvolve, regularizing with epsilon
 x= real(ifft(X)); % transform back to original domain
 x=x(1:1:Lx); % take valid portion of deconvolved vector
 out(counter,:)=x;
end

%averages out vector ignoring nans (same as nanmean function)
output = zeros(1, Lx);
for i = 1:Lx
    avgs = out(:,i);
    avgs = avgs(~isnan(avgs));
    output(1,i) = sum(avgs) / length(avgs);
end
output = output / max(output); % normalization
