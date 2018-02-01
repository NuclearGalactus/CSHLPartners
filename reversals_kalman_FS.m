%% experiment 1-  mimic of my reversal learning task
savepath = 'C:\Users\tcare\Desktop\Better Science';
load('TE.mat');
% initialize using Gershman code, prior variance for weights = 1,
% variance of process noise = 0.01, variance of measurement noise = 1
param = KTD_defparam;
rewardSize = .7;
punishSize = .3;
revFreq = 100; % block size (normal or reversed contingencies alternate)
sizer = size(TE.OdorValveIndex);
n = sizer(1); % number of trials
trials = (1:n)';
numCues = 2;
%cueType = TE.OdorValveIndex; % each trial contains either stimulus 1 or 2
X = zeros(n,numCues); % stimuli matrix nTrials x 2,  
X(TE.OdorValveIndex == 1, 1) = .8;
X(TE.OdorValveIndex == 2, 2) = .8;
datastuff = TE.csLicks.rate';
Odor1Trials = datastuff(TE.OdorValveIndex == 1);
Odor2Trials = datastuff(TE.OdorValveIndex == 2);


reversed = floor(trials / revFreq);
reversed = rem(reversed, 2) == 1;

    % randomly deliver reward on 80% of CS+ trials and punishment on 50% of CS-
    % trials
    rando = rand(n,1);
    r = zeros(n,1);
    %r(cueType == 1 & ~reversed & rando > 0.2) = rewardSize;
    %r(cueType == 1 & reversed & rando > 0.5) = punishSize;
   % r(cueType == 2 & reversed & rando > 0.2) = rewardSize;
   % r(cueType == 2 & ~reversed & rando > 0.5) = punishSize;
    rw = TE.ReinforcementOutcome;
    reward = zeros(size(rw));
    for i = (1:n)
       if(rw(i) == "Reward")
           reward(i) = 1;
       elseif(rw(i) == "Punish")
           reward(i) = 0.25;
       else
           reward(i) = 0.5;
       end
    end
    param.s = 2;
    param.q = 0.01;
    % run model
    model = kalmanRW(X,reward,param);
    % collect data
    rhat = zeros(n,1); % Predicted reward
    pe = zeros(n,1); % prediction error
    w = zeros(numCues,n); % weights
    Kn = zeros(numCues,n); % Kalman gain
    offDiag = zeros(n,1); % off diagonal term in posterior weight covariance matrix
    onDiag = zeros(2,n); % on diagonal terms
    for counter = 1:n
       rhat(counter) = model(counter).rhat;
       pe(counter) = model(counter).dt;
       w(:,counter) = model(counter).w0 * 2;
       Kn(:,counter) = model(counter).K;
       offDiag(counter) = model(counter).C(2,1); % covariance matrix is symmetric so bottom left or top right corner of 2,2 matrix are equivalent
       onDiag(:,counter) = [model(counter).C(1,1); model(counter).C(2,2)];
    end
    
    
    ensureFigure('Post and Pre Reversal Data');
    revTrials = find(TE.BlockChange)';
    subplot(2,2,1)
    numTrials = 100;
    colormap jet;
    cmap = colormap;
    hold on;
    for i = 1:length(revTrials)
        range = revTrials(i):(revTrials(i) + numTrials);
        dataset = TE.csLicks.rate(range);
        logical = TE.OdorValveIndex == 1;
        plot(find(logical(range)),kalman_noise_filter(dataset(logical(range))), 'Color', cmap(ceil(i/length(revTrials) * 64),:));
    end
    hold off;
    axis([0 numTrials 0 6]);
    subplot(2,2,2)
    hold on;
    for i = 1:length(revTrials)
        range = (revTrials(i) - numTrials):(revTrials(i));
        dataset = TE.csLicks.rate(range);
        logical = TE.OdorValveIndex == 1;
        plot(find(logical(range)),kalman_noise_filter(dataset(logical(range))), 'Color', cmap(ceil(i/length(revTrials) * 64),:));
    end
    hold off;
        axis([0 numTrials 0 6]);
    subplot(2,2,3)
    hold on;
    for i = 1:length(revTrials)
        range = (revTrials(i)):(revTrials(i) + numTrials);
        plot(w(1,range) * 2, 'Color', cmap(ceil(i/length(revTrials) * 64),:));
    end
    hold off;
        axis([0 numTrials 0 6]);
    subplot(2,2,4)
    hold on;
    for i = 1:length(revTrials)
        range = (revTrials(i) - numTrials):(revTrials(i));
        plot(w(1,range) * 2, 'Color', cmap(ceil(i/length(revTrials) * 64),:));
    end
    hold off;
        axis([0 numTrials 0 6]);
%{
    ensureFigure('Reversal_KalmanRW', 1);
   %  figure;
subplot(3,2,1);
plot(trials, w(1:2,:)');
ylim = get(gca, 'YLim');
revTrials = find(TE.BlockChange)';
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Weights (mean)');
xlabel('Trials');
subplot(3,2,2);
plot(trials, onDiag');
ylim = get(gca, 'YLim');
%line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Weights (variance)');
xlabel('Trials');

subplot(3,2,3);

plot(find(TE.OdorValveIndex == 1),kalman_noise_filter(TE.csLicks.rate(TE.OdorValveIndex == 1)), 'LineWidth', 1);
hold on;
plot(find(TE.OdorValveIndex == 2),kalman_noise_filter(TE.csLicks.rate(TE.OdorValveIndex == 2)), 'LineWidth', 1);
hold off;
ylim = get(gca, 'YLim');
line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Licking (rate)');
xlabel('Trials');
axis([0 n 0 ylim(2)]);
subplot(3,2,4);
plot(trials, abs(pe));
ylim = get(gca, 'YLim');
ylabel('Prediction error');
xlabel('Trials');
subplot(3,2,5);
s = scatter(w(1,TE.OdorValveIndex == 1),Odor1Trials);
hold on;
s.MarkerEdgeColor = 'b';
p = polyfit(w(1,TE.OdorValveIndex == 1 & ~isnan(TE.csLicks.rate)),Odor1Trials(~isnan(Odor1Trials)),1);
r = polyval(p, 0:.1:4.8);
plot(0:.1:4.8,r, 'LineWidth',1);
hold off;
ylim = get(gca, 'YLim');
ylabel('Licking (C1)');
xlabel('Weights (C1)');
axis([0 4.8 0 10]);
subplot(3,2,6);
s = scatter(w(2,TE.OdorValveIndex == 2 & ~isnan(TE.csLicks.rate)), Odor2Trials(~isnan(Odor2Trials)));
hold on;
s.MarkerEdgeColor = 'r';
p = polyfit(w(2,(TE.OdorValveIndex == 2 & ~isnan(TE.csLicks.rate))),Odor2Trials(~isnan(Odor2Trials)),1);
r = polyval(p, 0:.1:4.8);
plot(0:.1:4.8,r, 'LineWidth',1);
hold off;
ylim = get(gca, 'YLim');
%line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Licking (C2)');
xlabel('Weights (C2)');
axis([0 4.8 0 4]);
%{
ensureFigure("Brain Data");
subplot(2,2,1);
acetData = kalman_noise_filter(TE.phPeakMean_us(2).data);
plot(acetData(1,:));
hold on;
plot(acetData(2,:));
hold off;
ylim = get(gca, 'YLim');
%line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Dopamine');
xlabel('Trials');
axis([0 n -1 5]);
subplot(2,2,2);
acetData = kalman_noise_filter(TE.phPeakMean_us(1).data);
plot(acetData(1,:));
hold on;
plot(acetData(2,:));
hold off;
ylim = get(gca, 'YLim');
%line(repmat(revTrials, 2, 1), repmat(ylim', 1, length(revTrials)), 'Color', 'g');
ylabel('Acetylcholine');
xlabel('Trials');
axis([0 n -1 3]);
subplot(2,2,3);
s = scatter(w(1,:), acetData(1,:));
hold on;
s.MarkerEdgeColor = 'b';
p = polyfit(w(1,:),acetData(1,:),1);   % p returns 2 coefficients fitting r = a_1 * x + a_2
r = polyval(p, 0:.1:4.8);
plot(0:.1:4.8,r, 'LineWidth',1);
hold off;
ylim = get(gca, 'YLim');
ylabel('Licking (C1)');
xlabel('Weights (C1)');
axis([0 4.8 0 10]);
subplot(2,2,4);
weights = w(1);
s = scatter(w(1,:), acetData(1,:));
hold on;
s.MarkerEdgeColor = 'r';
p = polyfit(w(1,:),acetData(1,:),1);
r = polyval(p, 0:.1:4.8);
plot(0:.1:4.8,r, 'LineWidth',1);
hold off;
ylim = get(gca, 'YLim');
ylabel('Acetylcholine (C1)');
xlabel('Weights (C1)');
axis([0 4.8 0 4]);
%}
saveas(gcf, fullfile(savepath, 'Reversals_KalmanRW'), 'fig');
saveas(gcf, fullfile(savepath, 'Reversals_KalmanRW'), 'jpeg');
    %}