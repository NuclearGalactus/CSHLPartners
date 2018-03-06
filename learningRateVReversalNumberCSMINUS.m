%% learning rate v reversal number

learningRates = zeros(1,6);
path = 'C:\Users\tcare\Documents\Github\CSHLPartners\';
reversalNumbers = zeros(1,35);
load('AR.mat');
viewRange = -20:40;
reversalNumbers(1) = 1;
firstHalf = AR.csPlus;
secondHalf = AR.csMinus;
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


 figure;
 for i = 1:35
     subplot(7,5,i)
     hold on;
     plot(data(i,:));
     plot(smoothdata(data(i,:),'movmean',5), 'Color','r', 'LineWidth',2);
     %title("Reversal " + reversalNumbers(i));
     line([reversalPoint reversalPoint], [5 0], 'Color', [0 0 0], 'LineWidth', 1.5);
     line([tenPercentPoints(i) tenPercentPoints(i)], [5 0], 'Color', [0 0 0], 'LineWidth', 1.5);
     hold off;
 end

all = figure;
data = [firstHalf.csLicks.before secondHalf.csLicks.after];
dataRange = (viewRange) + size(firstHalf.csLicks.before,2);
for i = 1:6
    subplot(3,2,i);
    hold on;
    title("Reversal " + i);
    lr = getLearningRateCSMINUS(firstHalf,secondHalf,(reversalNumbers == i)',-4);
    meanData = nanmean(data(reversalNumbers == i,:));
    reversalPoint = size(firstHalf.csLicks.before,2);
    normalizer = mean(meanData(reversalPoint + (-20:-10)));
    plot(viewRange, nanmean(data(reversalNumbers == i,dataRange)) / normalizer, 'Color', 'g', 'LineWidth', 1.5);
    plot(viewRange, lr.graph(dataRange), 'Color', 'b', 'LineWidth',1.5); 
    reversals = find(reversalNumbers == i);
    for j = 1:length(reversals)
        plot(viewRange, data(reversals(j),dataRange) / normalizer, 'Color', 'r');
    end
    hold off;
    
    learningRates(1,i) = lr.value;
end
%saveas(all,path + "ReversalsCSPLUS.jpg");
%saveas(all,path + "ReversalsCSPLUS.fig");
learns = figure;
plot(learningRates, 'Color', 'r');
xticks([1 2 3 4 5 6]);
axis([1 6 0 1]);
title('Learning Rate over Multiple Reversals');
xlabel('Reversal Number');
ylabel('Optimal Learning Rate');
saveas(learns, path + "ReversalVLearningRateCSPLUS.fig");
saveas(learns, path + "ReversalVLearningRateCSPLUS.jpg");

