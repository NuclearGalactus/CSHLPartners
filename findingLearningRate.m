%% finding learning rate with correlation
load('AR.mat');
param = KTD_defparam;
rewardSize = .7;
punishSize = .3;
revFreq = 100; % block size (normal or reversed contingencies alternate) % number of trials
numCues = 2;


%X = zeros(n,numCues); % stimuli matrix nTrials x 2,  
%X(TE.OdorValveIndex == 1, 1) = .8;
%X(TE.OdorValveIndex == 2, 2) = .8;
mainDataPlus = [AR.csMinus.csLicks.before AR.csPlus.csLicks.after];
mainDataMinus = [AR.csPlus.csLicks.before AR.csMinus.csLicks.after];
n = size(mainDataPlus,1);
trials = (1:n)';
corrs = zeros(100,1);
biggestCorr = 0;
lastCorr = 0;


value = zeros(2,n);
totalCorr = 0;
rewards = [AR.csMinus.ReinforcementOutcome.before AR.csPlus.ReinforcementOutcome.after];
valves = [AR.csMinus.OdorValveIndex.before AR.csPlus.OdorValveIndex.after];
for lr = 1:2 - 1
    totalCorr = 0;
    for reversal = 1:n
        data = mainDataPlus(reversal,:);
        notnans = find(~isnan(data));
        data = data(notnans(1):notnans(end));
        numTrials = length(data); 
        r = zeros(numTrials,1);
        rw = rewards(reversal,:);
        rw = rw(notnans(1):notnans(end));
        reward = zeros(size(rw));
        for i = (1:numTrials)
         if(rw(i) == "Reward")
               reward(i) = 1;
         elseif(rw(i) == "Punish")
                reward(i) = .4;
         else
             reward(i) = 0.5;
         end
         X = zeros(2,numTrials);
         X(1,valves(reversal,:) == 1) = 1;
         X(2,valves(reversal,:) == 2) = 1;
         param.s = lr;
    
        param.q = 0.01;
        model = kalmanRW(X,reward,param);
        param.s = lr;
         
        param.q = 0.01;
  
    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(n,1); % prediction error
        w = zeros(numCues,n); % weights
        Kn = zeros(numCues,n); % Kalman gain
        offDiag = zeros(n,1); % off diagonal term in posterior weight covariance matrix
        onDiag = zeros(2,n); % on diagonal terms
        value = zeros(1,n);
        for counter = 1:n
           rhat(counter) = model(counter).rhat;
           pe(counter) = model(counter).dt;
           w(:,counter) = model(counter).w0 * 2;
           sum(w,2)
           Kn(:,counter) = model(counter).K;
           offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
           onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
           value(1,counter) = model(counter).w0 * X(:,counter);
           
        end
         
        totalCorr = totalCorr + corr(value,data)
    
    end
    
    
end
    
   
    
   
   
    
    
end







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
figure;
subplot(1,1,1);
%plot(corrs);
%}