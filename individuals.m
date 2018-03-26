%%
load('AR.mat');
normValues = mean(AR.csPlus.csLicks.before(:, end - 30:end), 2);
AR.csPlus.csLicks.before = bsxfun(@rdivide, AR.csPlus.csLicks.before, normValues);
AR.csPlus.csLicks.after = bsxfun(@rdivide, AR.csPlus.csLicks.after, normValues);
AR.csMinus.csLicks.before = bsxfun(@rdivide, AR.csMinus.csLicks.before, normValues);
AR.csMinus.csLicks.after = bsxfun(@rdivide, AR.csMinus.csLicks.after, normValues);
firstHalf = AR.csPlus;
secondHalf = AR.csMinus;
mainData = [firstHalf.csLicks.before secondHalf.csLicks.after];
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
reversalPoint = size(firstHalf.csLicks.before,2);
reversalNumbers = zeros(35,1);
filenames = AR.csPlus.filename.before(:,end);
reversalNumbers(1) = 1;
for i = 2:35
    if(strncmpi(filenames(i),(filenames(i-1)),5))
        reversalNumbers(i) = reversalNumbers(i-1) + 1;
    else
        reversalNumbers(i) = 1;
    end
end
tenPercentPoints = zeros(35,1);
for i = 1:35
    smoothed = smoothdata(mainData(i,:),'movmean',3);
    %highVal = nanmean(smoothed((-10:0) + reversalPoint));
    lessThans = find(smoothed((0:50) + reversalPoint) <= (0.1));
    if(isempty(lessThans))
        tenPercentPoints(i) = -1;
    else
    tenPercentPoints(i) = lessThans(1);
    end
end
lr = linspace(0,0.3,30);
viewBefore = 20;
viewAfter = 40;
viewRange = -viewBefore:viewAfter;
errorRange = (0:40) + reversalPoint;
graphRange = reversalPoint + viewRange;
bestmodels = zeros(35, size(mainData,2));
bestlrs = zeros(35,1);
besterrors = nan(35,1);
for reversal = 1:35
   if(tenPercentPoints(reversal) == -1)
        continue;
   end
   errors = zeros(1,length(lr));
   for lrc = 1:length(lr)
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
        param = KTD_defparam();
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
        errorRange = (0:tenPercentPoints(reversal)) + reversalPoint;
        modeldata = output(errorRange);
        datadata = mainData(reversal,errorRange);
        currentError = sum(abs(modeldata(~isnan(datadata)) - datadata(~isnan(datadata))));
        if(isnan(besterrors(reversal)) || currentError < besterrors(reversal))
            besterrors(reversal) = currentError;
            bestlrs(reversal) = lr(lrc);
            bestmodels(reversal,:) = output;
        end
   end
end


figure;
for i = 1:35
    if(tenPercentPoints(i) == -1)
        continue;
    end
    subplot(7,5,i);
    hold on;
     line([(tenPercentPoints(i)) (tenPercentPoints(i))], [-2 2], 'LineWidth', 1.8);
    plot(viewRange, mainData(i,reversalPoint + viewRange), 'r');
    %plot(viewRange, mainData(i,reversalPoint + viewRange), 'r');
    plot(viewRange, bestmodels(i,reversalPoint + viewRange), 'Color', 'b');

    hold off;
    axis([-viewBefore tenPercentPoints(i) + 4 -.5 2]);
end

figure;
hold on;
scatter(reversalNumbers, bestlrs);
