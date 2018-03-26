
%plot(nanmean(TE.Photometry.data.ZS,1));
% path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';
path = 'C:\Users\Adam\Documents\Repos\CSHLPartners\Deconvolution\';
cd(path);
load('TE.mat');
cued_highValue_Reward = filterTE(TE, 'trialType', 1);
cued_lowValue_Reward = filterTE(TE, 'trialType', 4);
uncuedPunishment = filterTE(TE, 'trialType', 8);
kernal = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [-2 2], 'FluorDataField', 'ZS');



data =  TE.Photometry.data.ZS;
TE.Photometry.data.Fourier = zeros(size(data));
h = kernal.Avg; 
% h = ones(1, 101); % try other kernals?
% h(50) = 1;
for counter = 1:size(data, 1)    
    y = data(counter,:)';
    Lx=length(y)-length(h)+1;   
    Lx2=pow2(nextpow2(Lx));    
    Y=fft(y, Lx2);		   
    H=fft(h, Lx2);		   
    X=(Y')./H;        		    
    x=real(ifft(X, Lx2));      
    x=x(1:1:Lx);               
    x=x/max(abs(x));
    TE.Photometry.data.Fourier(1:length(x)) = x;
end

fourier = ensureFigure('Fourier_test', 1);
phPlotAverageFromTE(TE, cued_highValue_Reward, 1, 'FluorDataField', 'ZS', 'window', [-4 3], 'linespec', 'k'); hold on;
phPlotAverageFromTE(TE, cued_highValue_Reward, 1, 'FluorDataField', 'Fourier', 'window', [-4 3], 'linespec', 'r'); 

% plot(x, 'r');

saveas(fourier, fullfile(path, 'deconvolution_fourier_fitz.fig'));
saveas(fourier, fullfile(path, 'deconvolution_fourier_fitz.jpg'));