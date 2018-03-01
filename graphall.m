reversalPoint = size(AR.csPlus.csLicks.before,2);
 figure;
 for i = 1:35
     subplot(7,5,i)
     hold on;
     plot(data(i,:));
     plot(smoothdata(data(i,:),'movmean',5), 'Color','r', 'LineWidth',2);
     title("Reversal " + reversalNumbers(i));
     line([reversalPoint reversalPoint], [5 0], 'Color', [0 0 0], 'LineWidth', 1.5);
     hold off;
 end
