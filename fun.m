x = linspace(0.414,0.4145,30);


y1 = 2 - (x.^2 + x.^2);
y2 = 2 .* x + 2 .* x;

figure;
hold on;
plot(x,y1);
plot(x, y2);