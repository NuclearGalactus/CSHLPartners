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
convolved = zeros(length(individual),1);
h = kernal.Avg;
for i = 2:(length(individual)-1)
    currentVal = 0;
    for j = 0:2
        currentVal = currentVal + (h(kernal.xData == (j)) * individual(i - j + 1));
    end
    convolved(i,1) = currentVal;
end
summation = figure;
plot(convolved);

saveas(summation, fullfile(path, 'deconvolution_summation.fig'));
saveas(summation, fullfile(path, 'deconvolution_summation.jpg'));