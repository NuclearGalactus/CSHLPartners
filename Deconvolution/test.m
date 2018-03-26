delta = zeros(1, 20 * 6); 
delta(20*3) = 1;   
kernel = (0:20) / 10; 
out = convolvSumm(delta, 0:20, kernel);
%out = convolvFourier(delta, kernel);
plot(out);