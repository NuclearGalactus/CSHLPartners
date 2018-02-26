%% learning rate v reversal number
learningRates = zeros(1,6);
reversalNumbers = zeros(1,35);
load('AR.mat');
viewRange = 0:30;
reversalNumbers(1) = 1;
filenames = AR.csPlus.filename.before(:,end);
for i = 2:35
    if(strncmpi(filenames(i),(filenames(i-1)),5))
        reversalNumbers(i) = reversalNumbers(i-1) + 1;
    else
        reversalNumbers(i) = 1;
    end
end
for i = 1:6
    lr = getLearningRate(AR.csMinus,AR.csPlus,(reversalNumbers == i)',viewRange);
    learningRates(1,i) = lr;
end
figure;
plot(learningRates, 'Color', 'r');
xticks([1 2 3 4 5 6]);
axis([1 6 0 1]);
title('Learning Rate over Multiple Reversals');
xlabel('Reversal Number');
ylabel('Optimal Learning Rate');
data = [AR.csPlus.csLicks.before AR.csMinus.csLicks.after];
dataRange = viewRange + size(AR.csPlus.csLicks.before,2);
figure;
hold on;
for i = 1:6
    subplot(3,2,i);
    hold on;
    title("Reversal " + i);
    plot(viewRange, nanmean(data(reversalNumbers == i,dataRange)), 'Color', 'y', 'LineWidth', 1.5);
    reversals = find(reversalNumbers == i);
    for j = 1:length(reversals)
        plot(viewRange, data(reversals(j),dataRange), 'Color', 'r');
    end
    hold off;
end
hold off;