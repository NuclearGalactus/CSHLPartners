load('TE.mat');
path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';
%plot(nanmean(TE.Photometry.data.ZS,1));
cued_highValue_Reward = filterTE(TE, 'trialType', 1);
cued_lowValue_Reward = filterTE(TE, 'trialType', 4);
uncuedPunishment = filterTE(TE, 'trialType', 8);
kernal = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [0 2], 'FluorDataField', 'ZS');
figure;
data =  TE.Photometry.data.ZS(cued_highValue_Reward,:);
individual = data(1,:);
plot(individual);
convolved = convolvSumm(data,kernal.xData,kernal.Avg);
figure; plot(convolved);

saveas(summation, fullfile(path, 'deconvolution_summation.fig'));
saveas(summation, fullfile(path, 'deconvolution_summation.jpg'));