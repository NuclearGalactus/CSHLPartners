load('TE.mat');
%plot(nanmean(TE.Photometry.data.ZS,1));
path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';

cued_highValue_Reward = filterTE(TE, 'trialType', 1);
cued_lowValue_Reward = filterTE(TE, 'trialType', 4);
uncuedPunishment = filterTE(TE, 'trialType', 8);
kernelavg = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [0 2], 'FluorDataField', 'ZS');
data =  TE.Photometry.data.ZS;
%h = kernelavg.Avg; 
h = zeros(3,1) + 1;
h(1,1) = 1;
h(2,1) = 1;
h = h';
%h = -(0:19) + 20;
n = 2;
L = size(data,2);

Lx=size(data,2)-length(h)+1;
Lx2=pow2(nextpow2(Lx));
W = (hann(size(data,2)));
y = zeros(n,L);
Y = zeros(n,Lx2);
H= fft(h, Lx2);
X = zeros(n,Lx2);
out = zeros(n,Lx);
for counter = 1:n
 y(counter,:) = (data(counter,:));    % Find smallest power of 2 that is > Lx
 Y(counter,:)=fft(y(counter,:) .* W1', Lx2);		   % Fast Fourier transform		   % Fast Fourier transform
 X(counter,:)=Y(counter,:).*H;        		   % 
 x=real(ifft(X(counter,:), Lx2));      % Inverse fast Fourier transform
 x=x(1:1:Lx);               % Take just the first N elements
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
figure;
subplot(3,2,1); plot(y(1,:) .* W');
title("Data");
subplot(3,2,2); plot((h));
title("Kernel");
subplot(3,2,3); plot(real(Y(1,:)));
title("Fourier Data");
%subplot(3,2,4); plot(real(ifft(Y(1,:))));
%title("Inverse Fourier");
subplot(3,2,4); plot(real((H)));
title("Fourier Kernel");
subplot(3,2,5); plot(real((X(1,:))));
title("Fourier Multiplied");
subplot(3,2,6); plot((out(1,:)));
title("Final Output");
%saveas(fourier, fullfile(path, 'deconvolution_fourier.fig'));
%saveas(fourier, fullfile(path, 'deconvolution_fourier.jpg'));


