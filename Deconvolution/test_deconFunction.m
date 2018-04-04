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


out = deconv_Fourier(data, h);
full = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
avged = nanmean(data);
plot(avged / max(avged));
plot(out);
legend("Before", "After");
axis([0 180 -2 2]);
title("Final Output");