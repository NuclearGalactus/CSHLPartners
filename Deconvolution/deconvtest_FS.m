
%plot(nanmean(TE.Photometry.data.ZS,1));
 path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';
%path = 'C:\Users\Adam\Documents\Repos\CSHLPartners\Deconvolution\';
cd(path);
load('TE.mat');
cued_highValue_Reward = filterTE(TE, 'trialType', 1);
cued_lowValue_Reward = filterTE(TE, 'trialType', 4);
uncuedPunishment = filterTE(TE, 'trialType', 8);
kernal = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [0 2], 'FluorDataField', 'ZS');



data =  TE.Photometry.data.ZS;
TE.Photometry.data.Fourier = zeros(size(data));
h = kernal.Avg;
%h = h(19:end - 19);
h(h<0) = 0;
h = smooth(h, 10);
h = h ./ sum(h);

%h = (20 - (0:19)) / sum(0:19); 
%h = 1;
% h = ones(1, 101); % try other kernals?
% h(50) = 1;
for counter = 1:size(data, 1)   
    uz1 = data(counter,:)';
    uz2 = h;
    L2 = numel(uz2);
    L = numel(uz1);
% 
% % Time vector
    T = [-(L/2+1):1:(L/2-2)];
% 
% % L-point symmetric Hann window in the column vector W
    W = hann(L); 
    W2 = hann(L2);
% 
% % Multiply input signals, uz1 and uz2, with Hann window and take FFT in
% % order to make sure that the ends of signal in frequency domain match up
% % while keeping everything reasonably smooth; this greatly diminish
% % spectral leakage
    uz1f = fft(W.*uz1,L); 
    uz2f = fft(W2.*uz2,L);
    %uz1f = uz1f(20:end-20);
    %uz2f = uz2f(20:end-20);
    L = numel(uz1f);
    Stmp = real(ifft((uz1f.*conj(uz2f))./(uz2f.*conj(uz2f)+eps*mean(uz2f.*conj(uz2f)))));
    S = [Stmp((L/2):L); Stmp(1:(L/2)-1)];
    TE.Photometry.data.Fourier(counter, 1:length(S)) = S;
end

fourier = ensureFigure('Fourier_test', 1);
phPlotAverageFromTE(TE, cued_highValue_Reward, 1, 'FluorDataField', 'ZS', 'window', [-4 3], 'linespec', {'k'}); hold on;
phPlotAverageFromTE(TE, cued_highValue_Reward, 1, 'FluorDataField', 'Fourier', 'window', [-4 3], 'linespec', {'r'}); 

% plot(x, 'r');

saveas(fourier, fullfile(path, 'deconvolution_fourier_fitz.fig'));
saveas(fourier, fullfile(path, 'deconvolution_fourier_fitz.jpg'));