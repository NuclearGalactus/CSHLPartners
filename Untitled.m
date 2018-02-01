%%

ensureFigure("Kalman Test");
hold on;
colormap jet;
cmap = colormap;
for i = 1:100
plot(kalman_noise_filter(TE.csLicks.rate,i), 'Color', cmap(ceil(i/100 * 64), :));
end
hold off;