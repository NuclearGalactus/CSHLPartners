function bestlr = getLearningRateCSMINUS(firstHalf, secondHalf, reversals, pointBefore)
lr = linspace(0.01,.4,30);
corrs = zeros(1,length(lr));
param = KTD_defparam;
%number of reversals
reversalPoint = size(firstHalf.csLicks.before,2);
reversals = reversals == 1;
mainData = [firstHalf.csLicks.before(reversals,:) secondHalf.csLicks.after(reversals,:)];
tenPercentPoints = zeros(35,1);
for i = 1:size(find(reversals))
    smoothed = smoothdata(mainData(i,:),'movmean',3);
    highVal = nanmean(smoothed((-10:0) + reversalPoint));
    lessThans = find(smoothed((0:50) + reversalPoint) <= (0.02 * highVal));
    if(isempty(lessThans))
        tenPercentPoints(i) = -1;
    else
    tenPercentPoints(i) = lessThans(1) + reversalPoint;
   end
end
numTrials = size(mainData,2);
n = size(mainData,1);
rewards = [firstHalf.ReinforcementOutcome.before(reversals,:) secondHalf.ReinforcementOutcome.after(reversals,:)];
valves = [firstHalf.OdorValveIndex.before(reversals,:) secondHalf.OdorValveIndex.after(reversals,:)];
corrPerReversal = zeros(n,1);
for lrc = 1:length(lr)
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
             reward(i) = -.1;
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
        model = kalmanRW(X,reward,param,1);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(2,1); % prediction error
        w = zeros(2,numTrials); % weights
        Kn = zeros(2,numTrials); % Kalman gain
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
           output(1,counter) = (w(plus,counter));
        end
        errorReversals(reversal,valves(reversal,:) == plus) = abs(output(valves(reversal,:) == plus) - mainData(reversal,valves(reversal,:) == plus));
        models(reversal,:) = max(output,0);
        dataRange =(reversalPoint):tenPercentPoints(reversal);
        tempnans = find(~isnan(mainData(reversal,:))); 
        thisData = mainData(reversal,dataRange) / mean(mainData(reversal,tempnans(1:10)));
        thisOutput = output(dataRange);
        if(tenPercentPoints(reversal) == -1)
            corrPerReversal(reversal) = -1;
        else
            corrPerReversal(reversal) = sum(abs(thisOutput - (models(reversal,dataRange))));
            %corrPerReversal(reversal) = sum(abs((thisData(~isnan(thisData)) - thisOutput(~isnan(thisData)))));
        end
    end
    average = nanmean(mainData(:,dataRange));
    average(isnan(average)) = 0;
    corrs(lrc) = sum(abs(average - nanmean(models(:,dataRange))));
    
end
bestlr.value = lr(find(corrs == min(corrs), 1, 'first'));

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
             reward(i) = -.1;
          else
             %X(i,1) = 0;
             %X(i,2) = 0;
             reward(i) = 0;
          end
        end
        param.s = lr(1,lrc);
        param.q = 0.01;
        param.std = 1;
        param.lr = bestlr.value;
        if(isempty(find(X(:,2))) == 0)
            plus =2;
        else
            plus = 1;
        end
        model = kalmanRW(X,reward,param,1);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(2,1); % prediction error
        w = zeros(2,numTrials); % weights
        Kn = zeros(2,numTrials); % Kalman gain
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
bestlr.graph = nanmean(max(models,0));
%figure;
%plot(lr,corrs);
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
             reward(i) =-.1;
          else
             %X(i,1) = 0;
             %X(i,2) = 0;
             reward(i) = 0;
          end
        end
        param.s = lr(1,lrc);
        param.q = 0.01;
        param.std = 1;
        param.lr = .1;
        if(isempty(find(X(:,2))) == 0)
            plus =2;
        else
            plus = 1;
        end
        model = kalmanRW(X,reward,param,1);    
        rhat = zeros(n,1); % Predicted reward
        pe = zeros(2,1); % prediction error
        w = zeros(2,numTrials); % weights
        Kn = zeros(2,numTrials); % Kalman gain
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
        models(reversal,:) = max(output,0);
    end
bestlr.extra = max(output,0);