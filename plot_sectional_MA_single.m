clear all
open('D:\zhuomian\新代码\fig_sigma\infty.fig');
fig = gcf; 
ax = findobj(fig, 'Type', 'axes'); 
if numel(ax) > 1
    ax = ax(1); 
end
surf_objects = findobj(ax, 'Type', 'surface');

if isempty(surf_objects)
    error('No surface data');
else
    X_single = surf_objects(1).XData;
    Y_single = surf_objects(1).YData;
    Z_single = surf_objects(1).ZData; 
    X_avg = surf_objects(2).XData;
    Y_avg = surf_objects(2).YData;
    Z_avg = surf_objects(2).ZData; 
end
close(fig); 


%% figure: n = 100 section.
figure
[row, col] =find(X_avg==100);
plot(X_avg(col(1),:),Z_avg(col(1),:),'LineWidth',3,'Color',[0.20, 0.65, 0.70]) 
hold on
plot(X_single(col(1),:),Z_single(col(1),:),'LineWidth',3,'Color', [0.85, 0.55, 0.50])  
xline(100,'--','LineWidth',2)
xlabel('M','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Risk','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
lgd=legend('Model averaging','Single model','M=n' );
lgd.FontSize = 40;
lgd.FontName = 'Times New Roman';
legend('boxoff')
set(gca, 'LineWidth', 1, 'FontSize', 30, 'FontName', 'Times New Roman', 'Box', 'off')
xlim([0,400])
ylim([0,3])



%% figure: M = 100 section.
figure
[row, col] =find(X_avg==100);
plot(Y_avg(:,col(1)),Z_avg(:,col(1)),'LineWidth',3,'Color',[0.20, 0.65, 0.70]) 
hold on
plot(Y_single(:,col(1)),Z_single(:,col(1)),'LineWidth',3,'Color', [0.85, 0.55, 0.50])  
xline(100,'--','LineWidth',2)
xlabel('n','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Risk','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
lgd=legend('Model averaging','Single model' ,'n=M');
lgd.FontSize = 40;
lgd.FontName = 'Times New Roman';
legend('boxoff')
set(gca, 'LineWidth', 1, 'FontSize', 30, 'FontName', 'Times New Roman', 'Box', 'off')
xlim([0,400])
ylim([0,3])

figure
c=0:1
plot(Y_avg(:,col(1)),Z_avg(:,col(1)),'LineWidth',3,'Color',[0.20, 0.65, 0.70]) 
hold on
plot(Y_single(:,col(1)),Z_single(:,col(1)),'LineWidth',3,'Color', [0.85, 0.55, 0.50])