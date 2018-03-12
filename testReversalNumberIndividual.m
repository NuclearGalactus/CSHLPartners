%%
load('AR.mat');
AR.csPlus.csLicks.before = bsxfun(@rdivide, AR.csPlus.csLicks.before, normValues);
AR.csPlus.csLicks.after = bsxfun(@rdivide, AR.csPlus.csLicks.after, normValues);
AR.csMinus.csLicks.before = bsxfun(@rdivide, AR.csMinus.csLicks.before, normValues);
AR.csMinus.csLicks.after = bsxfun(@rdivide, AR.csMinus.csLicks.after, normValues);
% end kludge
firstHalf = AR.csPlus;
secondHalf = AR.csMinus;
viewBefore = 20;
viewAfter = 40;
viewRange = (-viewBefore):viewAfter;
filenames = AR.csPlus.filename.before(:,end);
data = [firstHalf.csLicks.before secondHalf.csLicks.after];
tenPercentPoints = zeros(35,1);
reversalPoint = size(firstHalf.csLicks.before,2);
for i = 2:35    
    if(strncmpi(filenames(i),(filenames(i-1)),5))
        reversalNumbers(i) = reversalNumbers(i-1) + 1;
    else
        reversalNumbers(i) = 1;
    end
end


%  figure;
%  for i = 1:35
%      subplot(7,5,i)
%      hold on;
%      plot(data(i,:));
%      plot(smoothdata(data(i,:),'movmean',5), 'Color','r', 'LineWidth',2);
%      %title("Reversal " + reversalNumbers(i));
%      line([reversalPoint reversalPoint], [5 0], 'Color', [0 0 0], 'LineWidth', 1.5);
%      line([tenPercentPoints(i) tenPercentPoints(i)], [5 0], 'Color', [0 0 0], 'LineWidth', 1.5);
%      hold off;
%  end

all = figure;
data = [firstHalf.csLicks.before secondHalf.csLicks.after];
dataRange = (viewRange) + size(firstHalf.csLicks.before,2);
for i = 1:6
    subplot(3,2,i);
    hold on;
    title("Reversal " + i);
    reversals = find(reversalNumbers == i);
    thisReversalsLearningRates = zeros(length(reversals),1);
    for reversal = 1:length(reversals)
        code = zeros(35,1);
        code(reversals(reversal)) = 1;
        lr = getLearningRateCSMINUS(firstHalf,secondHalf,code,0);
        plot(viewRange, data(reversals(reversal),dataRange), 'Color', 'r');
        thisReversalsLearningRates(reversal) = lr.value;
    end
    meanData = nanmean(data(reversalNumbers == i,:));
    reversalPoint = size(firstHalf.csLicks.before,2);
    plot(viewRange, nanmean(data(reversalNumbers == i,dataRange)), 'Color', 'g', 'LineWidth', 1.5);
    %plot(viewRange, lr.graph(dataRange), 'Color', 'b', 'LineWidth',1.5); 
    axis([-20 40 0 1]);
    %plot(viewRange, lr.extra(dataRange), 'Color', 'b', 'LineWidth',1.5);
    hold off;
    deviations(1,i) = std(thisReversalsLearningRates);
    learningRates(1,i) = mean(thisReversalsLearningRates);
end
saveas(all,fullfile(path, 'ReversalsCSMINUS.jpg'));
saveas(all,fullfile(path, 'ReversalsCSMINUS.fig'));
learns = figure;
plot(learningRates, 'Color', 'r');
xticks([1 2 3 4 5 6]);
axis([1 6 0 1]);
title('Learning Rate over Multiple Reversals');
xlabel('Reversal Number');
ylabel('Optimal Learning Rate');
saveas(learns, fullfile(path, 'ReversalVLearningRateCSMINUS.fig'));
saveas(learns, fullfile(path, 'ReversalVLearningRateCSMINUS.jpg'));

