%% finding learning rate with correlation
load('AR.mat');
param = KTD_defparam;
numCues = 2;
%ensureFigure('LearningRate',1);
figure;
%SET DATA SOURCES
firstHalf = AR.csPlus;
secondHalf = AR.csMinus;
mainData = [firstHalf.csLicks.before secondHalf.csLicks.after];
reversalPoint = size(firstHalf.csLicks.before, 2);
numTrials = size(mainData,2);
%number of reversals
n = size(mainData,1);
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
acetyl = [firstHalf.phPeakMean_cs_ch1.before secondHalf.phPeakMean_cs_ch1.after];
dopamine = [firstHalf.phPeakMean_cs_ch2.before secondHalf.phPeakMean_cs_ch2.after];
firstReversals = cellfun(@(x,y) ~strcmp(x(1:5), y(1:5)),   firstHalf.filename.after(1:end-1,1), firstHalf.filename.after(2:end,1));
firstReversals = [false; firstReversals];
%normalize data
for reversal = 1:n
    data = AR.csPlus.csLicks.before(reversal,:);
    notnans = find(~isnan(data) & data ~= 0);
    mainData(reversal,:)  = mainData(reversal,:) / (mean(data(notnans(end-9:end))));
end
%SET VIEW RANGES FOR GRAPH
viewBefore = 20;
viewAfter = 30;
viewRange = (1:(viewAfter + viewBefore + 1))-(viewBefore + 1);
%SET COLOR MAP
colormap jet;
cmap = colormap;
%SETUP LEARNING RATES
lr = linspace(0.15,0.25,10);
%lr = [.2];
errorReversals = zeros(n,numTrials);
errorGraphs = zeros(length(lr),length(viewRange));
corrs = zeros(length(lr),1);
%csPlus Odor
%[ 1 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]
for lrc = 1:length(lr)
    lrc
    datas = nan(n,numTrials);
    models = nan(n,numTrials);
    zscores = zeros(n,numTrials);
    totalCorr = 0;
    for reversal = 1:n
        data = mainData(reversal,:);
        notnans = find(~isnan(data));
        numTrials = length(data); 
        rw = rewards(reversal,:);
        reward = zeros(numTrials,1);
        X = zeros(numTrials,2);
        X(valves(reversal,:) == 1,1) = 1;
        X(valves(reversal,:) == 2,2) = 1;
        emptyCells = (cellfun(@isempty,rw));
        for i = (1:numTrials)
          if(emptyCells(i))
             reward(i) = 0;
          elseif(rw(i) == "Reward")
             reward(i) = 1;
          elseif(rw(i) == "Punish")
             reward(i) = -1;
          else
             %X(i,1) = 0;
             %X(i,2) = 0;
             reward(i) = 0;
          end
        end
        param.s = lr(1,lrc);
        param.q = 0.01;
        param.std = 1;
        param.lr = lr(1,lrc);
        if(isempty(find(X(:,2))) == 0)
            plus =2;
        else
            plus = 1;
        end
        model = kalmanRW(X,reward,param);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(numTrials,1); % prediction error
        w = zeros(numCues,numTrials); % weights
        Kn = zeros(numCues,numTrials); % Kalman gain
        offDiag = zeros(numTrials,1); % off diagonal term in posterior weight covariance matrix
        onDiag = zeros(2,numTrials); % on diagonal terms
        output = zeros(1,numTrials);
        for counter = 1:numTrials
           rhat(counter) = model(counter).rhat;
           pe(counter) = model(counter).dt;
           w(:,counter) = model(counter).w0;
           Kn(:,counter) = model(counter).K;
           offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
           onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
           output(1,counter) = w(plus,counter);
        end
        errorReversals(reversal,valves(reversal,:) == plus) = abs(output(valves(reversal,:) == plus) - mainData(reversal,valves(reversal,:) == plus));
        models(reversal,:) = output;
       
    end
    %avglast10 = nanmean(AR.csPlus.csLicks.before(~firstReversals,:));
    %avglast10 = nanmean(avglast10(reversalPoint-9:reversalPoint));
    meaned = nanmean(mainData(~firstReversals,:));
    normalized= meaned / avglast10;
    meanedModel = nanmean(models(~firstReversals,:));
    normalizedModel = meanedModel / nanmean(meanedModel(1,(1:10) + reversalPoint  + viewAfter - 10));
    plotModel = (meanedModel(viewRange + reversalPoint));
    plotData = (meaned(viewRange + reversalPoint));
    errorGraphs(lrc,:) =mean(errorReversals(reversal,:));
    corrs(lrc) = corr(plotModel(viewBefore:end)', plotData(viewBefore:end)');
    %subplot(2,1,1);
    %hold on;
    %plot(viewRange, smoothdata((plotModel - plotData).^2,'movmean',5),'Color',cmap((ceil((lrc/length(lr)) * 64)),:));
    %hold off;
    subplot(2,2,2);
    hold on;
    plotModel(plotModel < 0) = 0;
    plot(viewRange,plotModel, 'Color',cmap((ceil((lrc/length(lr)) * 64)),:),'LineWidth',.8);
    hold off;
    axis([-viewBefore viewAfter -1.2 1.2]);
    subplot(2,2,4);
    hold on;
    %domain = find(valves(reversal,(0:30) + reversalPoint) == plus);
    %range =  nanmean(errorReversals(:,valves(reversal,(0:30) + reversalPoint) == plus));
    plot(viewRange((viewBefore + 1):end),(plotModel((viewBefore + 1):end) - plotData((viewBefore + 1):end)).^2,'Color',cmap((ceil((lrc/length(lr)) * 64)),:));
    xlabel('Trials Relative to Reversal') % x-axis label
    ylabel('Squared Error')
    hold off;
end
%subplot(2,1,1);
%plot(lr, corrs);
acetylData = (nanmean(acetyl(:,viewRange + reversalPoint)));
dopamData = (nanmean(dopamine(:,viewRange + reversalPoint)));
subplot(2,2,1);
hold on;
plot(viewRange,acetylData / 2, 'Color','m');
plot(viewRange,dopamData / 2, 'Color','b');
xlabel('Trials Relative to Reversal') % x-axis label
ylabel('Activity of Neurotransmitter cells') % y-axis label
hold off;
axis([-viewBefore viewAfter -1.1 1.1]);
subplot(2,2,2);
hold on;
%plot(-viewBefore:viewAfter, 1, 'Color','b');
for reversal = 1:n
    notnans = ~isnan(mainData(reversal,:));
    plot(find(notnans) - reversalPoint,smoothdata(mainData(reversal,notnans),'movmean',3),'Color','g', 'LineWidth', 0.5);

end
title('Responses to New CS+ Odor');
xlabel('Trials Relative to Reversal') % x-axis label
ylabel('Lick Rate or Model Weights (Normalized)') % y-axis label
hold off;
axis([-viewBefore viewAfter 0 2]);
subplot(2,2,3);
hold on;
plot(lr,corrs);
xlabel('Learning Rate') % x-axis label
ylabel('Correlation')
%{

ensureFigure('Learning Rate Graph',1);
subplot(2,1,1);
hold on;
for i = 1:n
    xvals = (1:size(mainDataPlus,2)) - size(AR.csMinus.csLicks.before,2);
    plot(xvals(~isnan(mainDataPlus(i,:))),kalman_noise_filter(mainDataPlus(i,~isnan(mainDataPlus(i,:)))'), 'Color', 'g');
end
boundedline((1:size(mainDataPlus,2)) - size(AR.csMinus.csLicks.before,2) , nanmean(mainDataPlus), nansem(mainDataPlus), 'alpha');
ylabel("Anticipation (CS+)");
xlabel("Trials (Relative to Reversal)");
hold off;
subplot(2,1,2);
hold on;
for i = 1:size(mainDataMinus,1)
    xvals = (1:size(mainDataMinus,2)) - size(AR.csPlus.csLicks.before,2);
    plot(xvals(~isnan(mainDataMinus(i,:))),kalman_noise_filter(mainDataMinus(i,~isnan(mainDataMinus(i,:)))'), 'Color', 'g');
end
boundedline((1:size(mainDataMinus,2)) - size(AR.csPlus.csLicks.before,2) , nanmean(mainDataMinus), nansem(mainDataMinus), 'alpha'); 
ylabel("Anticipation (CS-)");
xlabel("Trials (Relative to Reversal)");
hold off;


lastCorr;
ensureFigure('learning Rate');
plot(corrs);
%}