clear all
% Delete the models where cq=1, with equal weights.
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
        if Mmax(i)>=n
            kq(n)=[];
        end
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
                        RV_limit(iii,jjj)=sigma2*50;
                        RB_limit(iii,jjj)=theta'*theta*50;
                    end
                end
            end
        end
        RV_limit=RV_limit+RV_limit';
        RB_limit=RB_limit+RB_limit';
        for ttt=1:M
            kq=s(ttt);
            if kq/n<1
                RV_limit(ttt,ttt)=sigma2*kq/(n-kq);
                RB_limit(ttt,ttt)=theta(kq+1:p)'*theta(kq+1:p)*n/(n-kq);
            elseif kq/n>1
                RV_limit(ttt,ttt)=sigma2/(kq/n-1);
                RB_limit(ttt,ttt)=theta(1:kq)'*theta(1:kq)*(kq-n)/kq+theta(kq+1:p)'*theta(kq+1:p)*kq/(kq-n);
            else
                RV_limit(ttt,ttt)=sigma2*50;
                RB_limit(ttt,ttt)=theta'*theta*50;
            end
        end
        for t=1:100
            %% random weights (uniform)
            % omega=rand(M,1);
            % omega=omega./sum(omega);
            %% random weights (Lognormal)
            r_lognorm = random('Lognormal', 0, 3, [M, 1]); % mu=0, sigma=1
            omega = r_lognorm / sum(r_lognorm);
            %% random weights (Sparse Dirichlet via Gamma)
            % r_sparse = random('Gamma', 0.05, 1, [M, 1]);
            % omega = r_sparse / sum(r_sparse);
            %% random weights (dominant)
            % dominant_ratio = 0.95;
            % num_dominant = 1;
            % omega = zeros(M, 1);
            % idx_dominant = randperm(M, num_dominant);
            % omega(idx_dominant) = dominant_ratio / num_dominant;
            % remaining_weight = (1 - dominant_ratio) / (M - num_dominant);
            % omega(omega == 0) = remaining_weight;
            % omega = omega / sum(omega);
            %% equal weights
            % omega=ones(M,1)./M;
            riskB(t)=omega'*RB_limit*omega;
            riskV(t)=omega'*RV_limit*omega;
            risk2(t)=omega'*(RB_limit+RV_limit)*omega;
        end
        RiskB(i,j)=mean(riskB);
        RiskV(i,j)=mean(riskV);
        Risk2(i,j)=mean(risk2);
    end
     fprintf('Completion progress:: %d/%d\n', j, i);
end

%% Figure: Out-of-sample bias and variance varies with sample size n and model number M (SNR=1)
% n_values = NN;
% M_values = Mmax;
% [M_grid, n_grid] = meshgrid(M_values, n_values);
% figure('Position', [100, 100, 800, 600], 'Color', 'white');
% % color1 = [0, 0.4470, 0.7410];  % blue
% % color2 = [0.8500, 0.3250, 0.0980];  % orange
% color1 = [0.55, 0.25, 0.45];  % Purple
% color2 = [0.45, 0.80, 0.55];  % Green
% surf1 = mesh(M_grid, n_grid, RiskB', ...
%     'FaceAlpha', 0.6, ...
%     'EdgeAlpha', 0.3, ...
%     'LineWidth', 0.8, ...
%     'FaceColor', color1, ...
%     'EdgeColor', color1 * 0.7, ...
%     'DisplayName', 'Predictive Bias', ...
%     'FaceLighting','gouraud');   % Out-of-sample bias
% hold on
% surf2 = mesh(M_grid, n_grid, RiskV', ...
%     'FaceAlpha', 0.6, ...
%     'EdgeAlpha', 0.3, ...
%     'LineWidth', 0.8, ...
%     'FaceColor', color2, ...
%     'EdgeColor', color2 * 0.7, ...
%     'DisplayName', 'Predictive Variance', ...
%     'FaceLighting','gouraud');    % Out-of-sample variance
% axis tight;
% grid on;
% set(gca, ...
%     'LineWidth', 1, ...
%     'FontSize', 30, ...
%     'FontName', 'Times New Roman', ...
%     'XColor', [0.2 0.2 0.2], ...
%     'YColor', [0.2 0.2 0.2], ...
%     'ZColor', [0.2 0.2 0.2], ...
%     'GridColor', [0.4 0.4 0.4], ...
%     'GridAlpha', 0.3, ...
%     'Box', 'on', ...
%     'Projection', 'perspective');
% box on;
% xlabel('M ','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
% ylabel('n ','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
% zlabel('Risk','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
% legend([surf1, surf2], 'Out-of-sample bias', 'Out-of-sample variance', ...
%     'Location', 'northeast', ...
%     'FontSize', 22, ...
%     'FontName', 'Times New Roman', ...
%     'Box', 'on');
% zlim([0,3])





%% Figure: Out-of-sample risk varies with sample size n and model number M (SNR=1)
n_values = NN;
M_values = Mmax;
R =  Risk2';
[M_grid, n_grid] = meshgrid(M_values, n_values);
figure('Position', [100, 100, 800, 600], 'Color', 'white');
mesh(M_grid, n_grid, R, ...
    'FaceAlpha', 0.7, ...
    'LineWidth', 0.5, ...
    'FaceColor', 'interp','FaceLighting','gouraud');
camproj perspective
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
xlabel('M ','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
ylabel('n ','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
zlabel('Risk','FontSize', 30,'LineWidth',5, 'FontName', 'Times New Roman');
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
caxis([0, 5]);
colorbar('Location', 'eastoutside', 'FontName', 'Times New Roman', 'FontSize', 16);
set(gcf, 'Renderer', 'OpenGL');
set(gca, 'SortMethod', 'depth');
lighting gouraud;
material dull;
zlim([0,5])