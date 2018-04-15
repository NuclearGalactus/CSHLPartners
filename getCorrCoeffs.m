function r = getCorrCoeffs(data, valves, rewards, lr)


r = zeros(size(data,1),1);
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
        data = input(lrc,:);
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
       
        w = zeros(2,numTrials); % weights
       
        output = zeros(1,numTrials);
        for counter = 1:numTrials
          
           w(:,counter) = model(counter).w0;
           
           output(1,counter) = (w(plus,counter));
        end
       
        dataRange =(reversalPoint):(reversalPoint + 40);
        tempnans = find(~isnan(mainData(reversal,:))); 
        thisData = mainData(reversal,dataRange) / mean(mainData(reversal,tempnans(1:10)));
        thisOutput = output(dataRange);
        
            corrPerReversal(reversal) = corr(thisData(~isnan(thisData))',thisOutput(~isnan(thisData))');
            %corrPerReversal(reversal) = sum(abs((thisData(~isnan(thisData)) - thisOutput(~isnan(thisData)))));
        r(lrc,1) = corr(output, input(lrc,1));
    
end

