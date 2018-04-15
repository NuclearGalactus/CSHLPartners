load('NewData/TE.mat');
normValues = mean(AR.csPlus.csLicks.before(:, end - 30:end), 2);
AR.csPlus.csLicks.before = bsxfun(@rdivide, AR.csPlus.csLicks.before, normValues);
AR.csPlus.csLicks.after = bsxfun(@rdivide, AR.csPlus.csLicks.after, normValues);
AR.csMinus.csLicks.before = bsxfun(@rdivide, AR.csMinus.csLicks.before, normValues);
AR.csMinus.csLicks.after = bsxfun(@rdivide, AR.csMinus.csLicks.after, normValues);
firstHalf = AR.csPlus;
secondHalf = AR.csMinus;
mainData = [firstHalf.csLicks.before secondHalf.csLicks.after];
mainData = mainData(:,5:(end-5));
rewards = [firstHalf.ReinforcementOutcome.before secondHalf.ReinforcementOutcome.after];
valves = [firstHalf.OdorValveIndex.before secondHalf.OdorValveIndex.after];
valves = valves(:,5:(end-5));
rewards = rewards(:,5:(end-5));
rl = size(rewards, 2);
rh = size(rewards, 1);
reward = strcmp(rewards, "Reward");
mainData(isnan(mainData)) = 0;
mainData = mainData(1,:);
valves = valves(1,:);
reward = reward(1,:);
licking = mainData(1,(valves == 1 & ~isnan(mainData)));
rewarding = reward(1,valves == 1 & ~isnan(mainData));
%B = polyfit(1:length(reward), mainData(1,:),7);r
%reward(1,1:(end-1))
figure;
hold on;

x = linspace(0,1,100);
actual = licking(6:end)';
input = [ones(length(rewarding) - 5,1) rewarding(6:end)' rewarding(5:(end-1))' rewarding(4:(end-2))' rewarding(3:(end-3))' rewarding(2:(end-4))' rewarding(1:(end-5))'];
b = regress(licking(6:end)',input);
calc = input * b;
plot(calc);
%r2 =  1 - sum((actual - calc).^2)/sum((actual - mean(actual)).^2);
%log = mnrfit(1:(length(reward) - 1),changes + 1);
% 
% y = zeros(length(x),1);
% for i = 1:length(r)
%     y = y + ((r(end - (i-1))) * ((x) .^ (i - 1)));
% end


