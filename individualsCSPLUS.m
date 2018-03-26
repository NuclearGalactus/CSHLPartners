%%
load('AR.mat');
path = 'C:\Users\tcare\Documents\Github\CSHLPartners\Final_Figures\';
normValues = mean(AR.csPlus.csLicks.before(:, end - 30:end), 2);
AR.csPlus.csLicks.before = bsxfun(@rdivide, AR.csPlus.csLicks.before, normValues);
AR.csPlus.csLicks.after = bsxfun(@rdivide, AR.csPlus.csLicks.after, normValues);
AR.csMinus.csLicks.before = bsxfun(@rdivide, AR.csMinus.csLicks.before, normValues);
AR.csMinus.csLicks.after = bsxfun(@rdivide, AR.csMinus.csLicks.after, normValues);
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
tenPercentPoints = zeros(35,1) + 20;
firstHalf = AR.csMinus;
secondHalf = AR.csPlus;
mainData = [firstHalf.phPeakMean_cs_ch2.before secondHalf.phPeakMean_cs_ch2.after];
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
reversalPoint = size(firstHalf.csLicks.before,2);
lr = linspace(0,1,30);
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
        
        model = kalmanRW(X,reward,param,0);    
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


all = figure;
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
saveas(all, fullfile(path, 'allreversalscsplus.fig'));
saveas(all, fullfile(path, 'allreversalscsplus.jpg'));
scattered = figure;
hold on;
scatter(reversalNumbers, bestlrs);
myfit = fit(reversalNumbers, bestlrs, 'poly1');
plot(myfit, 'predfun');
axis([1 6 0 1]);
saveas(scattered, fullfile(path, 'scattercsplus.fig'));
saveas(scattered, fullfile(path, 'scattercsplus.jpg'));

meanes = figure;
hold on;
p = polyfit(reversalNumbers, bestlrs, 1);
prediction = (p(1) * (1:6)) + p(2);
deviations = zeros(6,1);
means = zeros(6,1);
for i = 1:6
    deviations(i) = std(bestlrs(reversalNumbers == i));
    means(i) = mean(bestlrs(reversalNumbers == i));
end
errorbar(means, deviations);
hold off;
saveas(meanes, fullfile(path, 'meanscsplus.fig'));
saveas(meanes, fullfile(path, 'meanscsplus.jpg'));