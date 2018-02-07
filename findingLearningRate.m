%% finding learning rate with correlation
load('AR.mat');
param = KTD_defparam;
rewardSize = .7;
punishSize = .3;
revFreq = 100; % block size (normal or reversed contingencies alternate) % number of trials
numCues = 2;
figure;

%X = zeros(n,numCues); % stimuli matrix nTrials x 2,  
%X(TE.OdorValveIndex == 1, 1) = .8;
%X(TE.OdorValveIndex == 2, 2) = .8;
firstHalf = AR.csMinus;
secondHalf = AR.csPlus;
mainData = [firstHalf.csLicks.before secondHalf.csLicks.after];
reversalPoint = size(secondHalf.csLicks.before, 2);
numTrials = size(mainData,2);
n = size(mainData,1);
trials = (1:n)';

biggestCorr = 0;
lastCorr = 0;
yeetcounter = 0;
num = 2;
value = zeros(2,n);

totalCorr = 0;
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
viewBefore = 20;
viewAfter = 40;
viewRange = (1:(viewAfter + viewBefore))-(viewBefore);
corrdata = zeros(num,n);
%ensureFigure('LearningRate',1);
% for lr = [1 100]
colormap jet;
cmap = colormap;
lr = linspace(0.05,1,19);
%lr = [0.1 1];
zdata = zeros(n,100);
corrs = zeros(length(lr),1);
for lrc = 1:length(lr)
    datas = nan(n,numTrials);
    models = nan(n,numTrials);
    zscores = zeros(n,numTrials);
    totalCorr = 0;
    for reversal = 1:n
        data = mainData(reversal,:);
        notnans = find(~isnan(data));
        data = data(notnans(1):notnans(end));
        numTrials = length(data); 
        rw = rewards(reversal,:);
        rw = rw(notnans(1):notnans(end));
        reward = zeros(numTrials,1);
        for i = (1:numTrials)
          if(rw(i) == "Reward")
               reward(i) = 1;
          elseif(rw(i) == "Punish")
                reward(i) = -1;
          else
             reward(i) = 0;
          end
        end
        X = zeros(numTrials,2);
        X(valves(reversal,notnans(1):notnans(end)) == 1,1) = 1;
        X(valves(reversal,notnans(1):notnans(end)) == 2,2) = 1;
    
        param.q = 0.01;

        param.s = lr(1,lrc);
         
        param.q = 0.01;
        param.std = 1;
        param.lr = lr(1,lrc);
        model = kalmanRW(X,reward,param);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(numTrials,1); % prediction error
        w = zeros(numCues,numTrials); % weights
        Kn = zeros(numCues,numTrials); % Kalman gain
        offDiag = zeros(numTrials,1); % off diagonal term in posterior weight covariance matrix
        onDiag = zeros(2,numTrials); % on diagonal terms
        value = zeros(1,numTrials);
        for counter = 1:numTrials
           rhat(counter) = model(counter).rhat;
           pe(counter) = model(counter).dt;
           w(:,counter) = model(counter).w0;
           Kn(:,counter) = model(counter).K;
           offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
           onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
           value(1,counter) = X(counter,:) * (model(counter).w0);   
        end
%        datas(reversal,:) = mainData(reversal,:);
        models(reversal,notnans(1):notnans(end)) = value;
        
        xs = (1:numTrials) - reversalPoint;
        z = (zscore(data(~isnan(data) & data ~= 0)) - zscore(value(~isnan(data) & data ~= 0))).^2;
        %ztemp(reversal, ~isnan(datas(reversal,:)) & datas(reversal,:) ~= 0) = (nanzscore(datas(reversal,~isnan(datas(reversal,:)) & datas(reversal,:) ~= 0)) - nanzscore(models(reversal,~isnan(datas(reversal,:)) & datas(reversal,:) ~= 0))).^2;
        %zscores(reversal,:) = ztemp;
        if(reversal == 1)
           % plot(find(~isnan(data) & data ~= 0),zscore(data(~isnan(data) & data ~= 0)),'Color','r');
           % plot(find(~isnan(data) & data ~= 0),zscore(value(~isnan(data) & data ~= 0)));
           % plot(find(~isnan(data) & data ~= 0),z, 'Color', 'p');
        end
        %plot(find(~isnan(data) & data ~= 0),z);
        %if(lrc == 1)
         %   plot((zscore(data(~isnan(data) & data ~= 0)) - zscore(value(~isnan(data) & data ~= 0))).^2, 'Color', 'g');
       % elseif(lrc == 2)
        %    plot((zscore(data(~isnan(data) & data ~= 0)) - zscore(value(~isnan(data) & data ~= 0))).^2, 'Color', 'r');
       % end
        %currentdiff  = corr(value(~isnan(data) & data ~= 0)',data(~isnan(data) & data ~= 0)');
        %corrdata(lrc,reversal) = currentCorr;
       % totalCorr = totalCorr + currentCorr;
       
       
       
    end
    corrs(lrc) = corr(zscore(nanmean((models(:,viewRange + reversalPoint))))', zscore(nanmean(mainData(:,viewRange + reversalPoint)))');
    
    subplot(2,1,2);
    hold on;
    plot(viewRange,zscore(nanmean((models(:,viewRange + reversalPoint)))), 'Color',cmap((ceil((lrc/length(lr)) * 64)),:),'LineWidth',.8);
    
    hold off;
    axis([-viewBefore viewAfter -2 2]);
end
subplot(2,1,1);
plot(lr, corrs);

subplot(2,1,2);
hold on;
plot(viewRange,zscore(nanmean(mainData(:,viewRange + reversalPoint))), 'Color','g','LineWidth',.8);
hold off;
%plot(lr,mean(zdata));




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