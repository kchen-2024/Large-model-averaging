clear all
% double descent
% Nested model averaging
NN=2:2:400;
Mmax=2:2:400;
p=400;
sigma2=1;
snr=1;
%% Polynomial decay
alpha=0.1;
theta = ((1:p).^(-alpha-0.5)) * sqrt(2*alpha);
theta=theta';
theta=theta./norm(theta)*sqrt(sigma2*snr);
for j=1:length(NN)
    n=NN(j);
    for i=1:length(Mmax)
        kq=1:1:Mmax(i);
        % if Mmax(i)>=n
        %     kq(n)=[];
        % end
        s=kq;
        %% Calculate the limit of out-of-sample risk
        M=length(s);
        RB_limit=zeros(M,M);
        RV_limit=zeros(M,M);
        for iii=1:M
            for jjj=1:M
                if jjj>iii
                    ki=s(iii);
                    kj=s(jjj);
                    if ki<n && kj<n
                        RV_limit(iii,jjj)=sigma2*ki/(n-ki);
                        RB_limit(iii,jjj)=theta(kj+1:p)'*theta(kj+1:p)*n/(n-ki);
                    elseif ki<n && kj>n
                        RV_limit(iii,jjj)=sigma2*ki/(kj-ki);
                        RB_limit(iii,jjj)=theta(kj+1:p)'*theta(kj+1:p)*kj/(kj-ki)+theta(ki+1:kj)'*theta(ki+1:kj)*(kj-n)/(kj-ki);
                    elseif ki>n && kj>n
                        RV_limit(iii,jjj)=sigma2/(kj/n-1);
                        RB_limit(iii,jjj)=theta(1:ki)'*theta(1:ki)*(ki-n)/ki+theta(ki+1:kj)'*theta(ki+1:kj)+theta(kj+1:p)'*theta(kj+1:p)*kj/(kj-n);
                    else
                        RV_limit(iii,jjj)=sigma2*500;
                        RB_limit(iii,jjj)=theta'*theta*500;
                    end
                end
            end
        end
        RV_limit=RV_limit+RV_limit';
        RB_limit=RB_limit+RB_limit';
        for ttt=1:M
            kq=s(ttt);
            c_q=kq/n;
            if kq/n<1
                RV_limit(ttt,ttt)=sigma2*kq/(n-kq);
                RB_limit(ttt,ttt)=theta(kq+1:p)'*theta(kq+1:p)*n/(n-kq);
            elseif kq/n>1
                RV_limit(ttt,ttt)=sigma2/(kq/n-1);
                RB_limit(ttt,ttt)=theta(1:kq)'*theta(1:kq)*(kq-n)/kq+theta(kq+1:p)'*theta(kq+1:p)*kq/(kq-n);
            else
                RV_limit(ttt,ttt)=sigma2*500;
                RB_limit(ttt,ttt)=theta'*theta*500;
            end
        end
        omega=ones(M,1)./M;
        Risk2(i,j)=omega'*(RB_limit+RV_limit)*omega;
        %% Calculate the risk of a single model kM
        c=max(s)/n;
        if c>1
            Rx(i,j)=norm(theta(1:max(s)))^2*(1-1/c)+sigma2/(c-1)+norm(theta(1+max(s):end))^2*c/(c-1);
        elseif c==1
            Rx(i,j)=snr*500;
        else
            Rx(i,j)=sigma2*c/(1-c)+norm(theta(1+max(s):end))^2/(1-c);
        end
        fprintf('Completion progress:: %d/%d\n',i, j);
    end
end


%% Combined Graph: Comparison of Individual Maximum Model with Model Averaging
n_values = NN;
M_values = Mmax;
[M_grid, n_grid] = meshgrid(M_values, n_values);
figure('Position', [100, 100, 950, 700], 'Color', 'white');
custom_map = [
    0.00  0.00  0.25
    0.00  0.20  0.70
    0.00  0.50  0.85
    0.00  0.80  0.70
    0.40  0.95  0.50
    0.60  0.95  0.40
    0.80  0.95  0.30
    0.90  0.90  0.20
    0.95  0.75  0.10
    0.95  0.55  0.00
    0.95  0.35  0.00
    0.925 0.175 0.00
    0.90  0.00  0.00
    ];
colormap(custom_map);
caxis([0, 1.8]);
color_single = [0.75, 0.25, 0.20];
surf1 = surf(M_grid, n_grid, Risk2', ...
    'FaceAlpha', 1, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none', ...
    'CData', Risk2');
hold on

edge_color = color_single * 0.8; 
surf2 = surf(M_grid, n_grid, Rx', ...
    'FaceAlpha', 0.25, ...              % 保持你原来的 0.25 透明度
    'FaceColor', color_single, ...      % 保持你原来的纯色，绝无渐变
    'EdgeColor', edge_color, ...        % 【核心修改】开启网格线，用来勾勒曲面的折叠和弯曲
    'EdgeAlpha', 0.35, ...              % 【核心修改】给网格线设置半透明，让它隐隐约约，不刺眼
    'LineWidth', 0.5, ...               % 线条细一点，更显精致
    'FaceLighting', 'none');

%camlight('headlight');
lighting gouraud;

axis tight;
grid on;
set(gca, ...
    'LineWidth', 1, ...
    'FontSize', 30, ...
    'FontName', 'Times New Roman', ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2], ...
    'ZColor', [0.2 0.2 0.2], ...
    'GridColor', [0.4 0.4 0.4], ...
    'GridAlpha', 0.3, ...
    'Box', 'on', ...
    'Projection', 'perspective');
box on;
zlim([0, 8]);
xlabel('M', 'FontSize', 26, 'FontName', 'Times New Roman');
ylabel('n', 'FontSize', 26, 'FontName', 'Times New Roman');
zlabel('Risk', 'FontSize', 26, 'FontName', 'Times New Roman');

caxis([min(Risk2(:)), 1.8]);
cb = colorbar('Location', 'eastoutside', 'FontName', 'Times New Roman', 'FontSize', 16);
cb.Label.String = 'Model averaging risk';
cb.Label.FontSize = 18;

lgd_p1 = patch([0 1 1 0], [0 0 1 1], [0 0 0 0], ...
    'FaceColor', custom_map(2,:), 'EdgeColor', 'none');
lgd_p2 = patch([0 1 1 0], [0 0 1 1], [0 0 0 0], ...
    'FaceColor', color_single, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

legend([lgd_p1, lgd_p2], {'Model averaging', 'Single model'}, ...
    'Location', 'northeast', ...
    'FontSize', 22, ...
    'FontName', 'Times New Roman', ...
    'Box', 'on');
zlim([0,3])