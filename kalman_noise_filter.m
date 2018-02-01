function betterdata = kalman_noise_filter(inputData, varargin)
%this is a github test
if length(varargin) < 1
    intensity = 1;
else
    intensity = varargin{1}; 
end

[T,N] = size(inputData);
if(T > N)
    inputData = inputData';
    [~,T] = size(inputData);
end
olddata = inputData;
olddata(isnan(olddata)) = 0;
betterdata = zeros(size(olddata));
variance = 0;
noise = intensity;
q = .01;
betterdata(1,1) = olddata(1,1);
K = 0;
for i = (2:T)
    estimated = betterdata(1, i-1);
    variance = variance + q;
    K = variance / (variance + noise);
    perr = (olddata(1, i) - estimated);
    betterdata(1,i) = estimated + (K * perr);
    variance = (1 - K) * variance;
end



%{
 ensureFigure('Kalman test', 1);
 subplot(1,1,1); 
 plot(olddata(1,:));
 hold on;
 axis([0 1000 0 10]);
 plot(betterdata(1,:));
 
 axis([0 1000 0 10])
 plot(smooth(olddata(1,:),11), 'g');
hold off;
%}