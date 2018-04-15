function bestlr = findBestLearningRateRW(mainData, cues, lr, rewards, ranges, logical, initial)
%reversals is a logical that will select which rows to include
reversals = logical;
param = KTD_defparam; %establish parameters for the learning algorithm
mainData = mainData(reversals,:); %establish actual dataset (must be analagous to predicted reward)
numTrials = size(mainData,2); %max length of reversals, used for padding
n = size(mainData,1); %number of reversals to iterate through
rewards = rewards(reversals,:); %defines how much reward is delivered during each trial
valves = cues(reverals,:); %defines which cues are presented during each trial
dataRange = ranges; %where to look in each 
bestlr.deviation = zeros(length(lr),1);
for lrc = 1:length(lr)
    models = nan(n,numTrials);
    errors = zeros(n,numTrials);
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
        model = kalmanRW(X,reward,param,0);    
        w = zeros(2,numTrials); % weights
        output = zeros(1,numTrials);
        for counter = 1:numTrials
           w(:,counter) = model(counter).w0;
           output(1,counter) = (w(plus,counter));
        end
        errorReversals(reversal,valves(reversal,:) == plus) = abs(output(valves(reversal,:) == plus) - mainData(reversal,valves(reversal,:) == plus));
        models(reversal,:) = max(output,0);
    end
    average = nanmean(mainData(:,dataRange));
    average(isnan(average)) = 0;
    errors(lrc) = sum(abs(average - nanmean(models(:,dataRange))));    
end
bestlr.value = lr(find(errors == min(errors), 1, 'first'));

%figure;
%plot(lr,corrs);
