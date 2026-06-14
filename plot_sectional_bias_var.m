clear all
open('D:\zhuomian\新代码\fig_sigma\bias_var_delete1.fig');
fig = gcf; 
ax = findobj(fig, 'Type', 'axes');
if numel(ax) > 1
    ax = ax(1); 
end
xlabel_text = get(get(ax, 'XLabel'), 'String');
ylabel_text = get(get(ax, 'YLabel'), 'String');
data_objects = findobj(ax, '-depth', 1, ...
    {'-property', 'XData', '-property', 'YData', '-property', 'ZData'});
if isempty(data_objects)
    error('No surface data');
else
    for i = 1:numel(data_objects)
        obj = data_objects(i);
    end
    obj = data_objects(1);
    X = obj.XData;
    Y = obj.YData;
    Z = obj.ZData;
    obj = data_objects(2);
    X2 = obj.XData;
    Y2 = obj.YData;
    Z2 = obj.ZData;
end
close(fig);



%% figure: n = 100 section.
[row, col] =find(X==100);
plot(X2(col(1),:),Z2(col(1),:),'LineWidth',3,'Color',[0.55, 0.25, 0.45]) 
hold on
plot(X(col(1),:),Z(col(1),:),'LineWidth',3,'Color', [0.45, 0.80, 0.55])  
xline(100,'--','LineWidth',2)
xlabel('M','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Risk','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
lgd=legend('Out-of-sample bias','Out-of-sample variance','M=n' );
lgd.FontSize = 40;
lgd.FontName = 'Times New Roman';
legend('boxoff')
set(gca, 'LineWidth', 1, 'FontSize', 30, 'FontName', 'Times New Roman', 'Box', 'off')
xlim([0,400])
ylim([0,3])



%% figure: M = 100 section. 
figure
[row, col] =find(X==100);
plot(Y2(:,col(1)),Z2(:,col(1)),'LineWidth',3,'Color',[0.55, 0.25, 0.45]) 
hold on
plot(Y(:,col(1)),Z(:,col(1)),'LineWidth',3,'Color', [0.45, 0.80, 0.55])  
xline(100,'--','LineWidth',2)
xlabel('n','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
ylabel('Risk','LineWidth',5, 'FontName', 'Times New Roman', 'FontSize', 10);
lgd=legend('Out-of-sample bias','Out-of-sample variance' ,'n=M');
lgd.FontSize = 40;
lgd.FontName = 'Times New Roman';
legend('boxoff')
set(gca, 'LineWidth', 1, 'FontSize', 30, 'FontName', 'Times New Roman', 'Box', 'off')
xlim([0,400])
ylim([0,3])
