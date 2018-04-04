load('TE.mat');
%plot(nanmean(TE.Photometry.data.ZS,1));
%path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';
path = 'C:\Users\tcare\OneDrive\Documents\GitHub\CSHLPartners\Deconvolution\';
cued_highValue_Reward = filterTE(TE, 'trialType', 1);
cued_lowValue_Reward = filterTE(TE, 'trialType', 4);
uncuedPunishment = filterTE(TE, 'trialType', 8);
kernel = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [0 2], 'FluorDataField', 'ZS');
data =  TE.Photometry.data.ZS;
h = kernel.Avg; 
%hx = linspace(0,2,20);
%h = exp(-1 * hx);
%h = zeros(3,1) + 1;
%h(1,1) = 1;
%h(2,1) = 1;
%h = h';
%h = -(0:19) + 20;
n = size(data,1);
L = size(data,2);
eps = 0.1;
Lx=size(data,2)-length(h)+1;
Lx2=pow2(nextpow2(Lx));
W = (blackman(size(data,2)));
y = zeros(n,L) + nanmean(data);
Y = zeros(n,Lx2);
H= fft(h, Lx2);
X = zeros(n,Lx2);
out = zeros(n,Lx);
for counter = 1:n
 y(counter,:) = data(counter,:);  % Find smallest power of 2 that is > Lx
 Y(counter,:)=fft(y(counter,:) .* W', Lx2);		   % Fast Fourier transform		   % Fast Fourier transform
 X(counter,:) = (Y(counter,:).*conj(H))./(H.*conj(H)+eps*mean(H.*conj(H))); 
 x= real(ifft(X(counter,:)));
 x=x(1:1:Lx);
 out(counter,:)=x;
end
% uz1 = data(1,:)';
% uz2 = kernelavg.Avg';
% eps = 0.1; 
% L = numel(uz1);
%  
% 
% T = [-(L/2+1):1:(L/2-2)];
% 
% 
% W = hann(L); 
% W2 = hann(numel(uz2));
% 
% uz1f = fft(W.*uz1,L); 
% uz2f = fft(W2.*uz2,L);
%  
% 
% %Stmp = real(ifft((uz1f.*conj(uz2f))));%./(uz2f.*conj(uz2f)+eps*mean(uz2f.*conj(uz2f)))));
% Stmp = real(ifft(uz1f ./ uz2f));
% S = [Stmp(L/2:L); Stmp(3:L/2-1)];
full = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,2,1); plot(nanmean(y));
title("Averaged Data");
subplot(3,2,2); plot((h));
title("Kernel");
subplot(3,2,3); plot(real(nanmean(Y)));
title("Averaged Data after Transformation");
%subplot(3,2,4); plot(real(ifft(Y(1,:))));
%title("Inverse Fourier");
subplot(3,2,4); plot(real((H)));
title("Kernel after Transformation");
subplot(3,2,5); plot(real(nanmean(X)));
title("Transforms Multiplied");
subplot(3,2,6); hold on;plot(nanmean((out)) / max(nanmean((out))));
plot(nanmean((y)) / max(nanmean(y)));
legend("Before", "After");
axis([0 180 -2 2]);
title("Final Output");
saveas(full, fullfile(path, 'deconvolution_fourier.fig'));
saveas(full, fullfile(path, 'deconvolution_fourier.jpg'));


