figure;
hold on;
load('TE.mat')
LNL_conditions;
licking = TE.csLicks.rate;
rewards = double(~rewardTrials(Odor1Trials,1));
licking = licking(Odor1Trials,1);


actual = licking(6:end);
log = ~isnan(actual);
actual = actual(log,1);
%ones(length(rewarding) - 5,1)
input = [ones(length(rewards) - 5,1) rewards(5:(end-1)) rewards(4:(end-2)) rewards(3:(end-3))];
input = input(log,:);
% rewarding(2:(end-4)) rewarding(1:(end-5))];
[b, bint, residuals, rint, stats] = regress(actual,input);
calc = (input * b);
plot(calc);
plot(actual);
r = corr(calc,actual)