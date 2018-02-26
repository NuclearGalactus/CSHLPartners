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
for i = 1:5
    lr = getLearningRate(AR.csPlus,AR.csMinus,(reversalNumbers == i)',viewRange);
    learningRates(1,i) = lr;
end
figure;
plot(learningRates, 'Color', 'r');
xticks([1 2 3 4 5]);
axis([1 5 0 0.2]);
title('Learning Rate over Multiple Reversals');
xlabel('Reversal Number');
ylabel('Optimal Learning Rate');