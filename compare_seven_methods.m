clear all
rng(7); 
n = 300;
p = 1000;
alpha = 0.5;
K = 500;
% M = round(3 * (n^(1/3))); 
M = round(0.9*n);
s = 1:M;
sigma_2 = 1;

%% Generate symmetric positive definite matrix
A1 = diag(rand(p,1));
A2 = orth(rand(p,p));
A = A2'*A1*A2;   

%% Set regression coefficients
R = 0.1:0.1:0.9;
rl = length(R);
c = sqrt( ( R*sigma_2 ) ./ (1 - R) );
the = ((1:p).^(-alpha-0.5)) * sqrt(2*alpha);

for r = 1:rl
    theta = c(r) * the(:); % Regression coefficient
    for g = 1:K

        %% Generating training data
        Xt = mvnrnd(zeros(1,p),  A , n);%   eye(p)
        Xt(:,1) = 1;
        e = randn(n,1) * sqrt(sigma_2);
        mu = Xt * theta;
        y = mu + e;
        X = Xt(:,1:M);
        %% Generating test data
        N_test=n;
        X_test = mvnrnd(zeros(1,p), A  , N_test);   %  eye(p)
        X_test(:,1) = 1;
        e_test=sqrt(sigma_2)*randn(N_test, 1);
        mu_test=X_test*theta;
        Y_test=mu_test+e_test;

        %% JMA
        ehat_jack = zeros(n,M);
        for k = 1:M
            X_sub = X(:,1:k);
            Pmm= X_sub*inv(X_sub'*X_sub)*X_sub';
            Dmm=diag(1./diag(eye(n)-Pmm));
            Pm_tilde=Dmm*(Pmm-eye(n))+eye(n);
            mu_tilde(:, k) =  Pm_tilde*y;
        end
        Y = repmat(y, 1, M);
        ehat_jack = Y - mu_tilde;
        % Quadratic programming solves the weights
        c1 = 2 * (ehat_jack' * ehat_jack) / n;
        if rank(c1) < M
            c1 = c1 + eye(M) * 1e-10;
        end
        c2 = zeros(M,1);
        Aeq = ones(1,M);
        beq = 1;
        lb = zeros(M,1);
        ub = ones(M,1);
        w0 = ones(M,1)./M;
        options = optimoptions('quadprog', 'Display', 'off');
        w_jack = quadprog(c1, c2, [], [], Aeq, beq, lb, ub, w0 ,options);
        w_jack(w_jack < 0) = 0;
        w_jack = w_jack / sum(w_jack);

        %% MMA
        % noise variance estimated by largest model
        % mm=M;
        % H = Xt(:,1:mm) *pinv (Xt(:,1:mm)'*Xt(:,1:mm)) * Xt(:,1:mm)';
        % e_large =(eye(n) - H) * y;
        % sigma = e_large' *e_large / (n - mm);
       
        thetahat = zeros(M,M);
        for k = 1:M
            Xk = Xt(:,1:s(k));
            thetahat(1:s(k), k) = Xk \ y;
        end
        muhat = X * thetahat;
        SSR = sum((y - muhat).^2, 1);
        Y = repmat(y, 1, M);
        ehat = Y - muhat;
        a1 = ehat' * ehat;
        if rank(a1) < M
            a1 = a1 + eye(M) * 1e-10;
        end
        % noise variance estimated by BIC-selected model
        bic = n * log(SSR)- n * log(n)+ s * log(n);
        [~, bicch] = min(bic);
        sigma = SSR(bicch) / (n - s(bicch));
        a2 = sigma * s(:);
        w = quadprog(a1, a2, [], [], Aeq, beq, lb, ub, w0, options);
        w = max(w, 0);
        w = w ./ sum(w);

       %% AIC/BIC
        aic = n * log(SSR) + s * 2 - n * log(n);
        [~, aicch] = min(aic);
        wa = exp(-0.5 * aic);
        wa = wa / sum(wa);
        wb = exp(-0.5 * bic);
        wb = wb / sum(wb);

        %% Large MA
        Rn1=zeros(M,1);
        An1=zeros(M,1);
        for t=1:M
            c_q=s(t)/n;
            if c_q<1
                Rn1(t)=sigma*c_q/(1-c_q);
                An1(t)=sigma*c_q;
            end
        end
        Rn2=zeros(M,M);
        An2=zeros(M,M);
        for iii=1:M
            for jjj=1:M
                if jjj>iii
                    ki=s(iii);
                    kj=s(jjj);
                    if ki<n && kj<n
                        Rn2(iii,jjj)=sigma*min(ki,kj)/(n-min(ki,kj));
                        An2(iii,jjj)=sigma*max(ki,kj)/n;
                    end
                end
            end
        end
        An=diag(An1)+An2+An2';
        Rn=diag(Rn1)+Rn2+Rn2';
        % Choice of the regularization parameter
        sm=max(s);
        delta=(n)/(n-sm);
        Q_V = delta*diag(Rn1)+Rn;
        if rank(Q_V) < size(ehat,2)
            Q_V = Q_V + eye(M)*1e-10;
        end
        % Solve the quadratic programming problem
        w_LaMA = quadprog(2*(a1+n*Q_V+n*An), zeros(M, 1), [], [], Aeq, beq, lb, ub, w0, options);
        w_LaMA = max(w_LaMA, 0);
        w_LaMA = w_LaMA ./ sum(w_LaMA);
        thetahatLaMA = thetahat * w_LaMA;

        %% Calculate the estimators of each method
        thetahataic = thetahat(:, aicch);
        thetahatbic = thetahat(:, bicch);
        thetahatMMA = thetahat * w;
        thetahatJMA = thetahat * w_jack;
        thetahatSAIC = thetahat * wa';
        thetahatSBIC = thetahat * wb';
        %% Baseline risk
        mumu = repmat(mu, 1, M);
        Risk0(g,:) = sum((mumu - muhat).^2, 1);
        muhat_test = X_test(:,1:M) * thetahat;
        mumu_test = repmat(mu_test, 1, M);
        Risk0_test(g,:) = sum((mumu_test  - muhat_test ).^2, 1);
        %% training risk
        Riskaic(g,1) = sum((mu - X * thetahataic).^2);
        Riskbic(g,1) = sum((mu - X * thetahatbic).^2);
        RiskMMA(g,1) = sum((mu - X * thetahatMMA).^2);
        RiskJMA(g,1) = sum((mu - X * thetahatJMA).^2);
        RiskSAIC(g,1) = sum((mu - X * thetahatSAIC).^2);
        RiskSBIC(g,1) = sum((mu - X * thetahatSBIC).^2);
        RiskLaMA(g,1) = sum((mu - X*thetahatLaMA).^2);
        %% test risk
        Riskaic_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahataic).^2);
        Riskbic_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahatbic).^2);
        RiskMMA_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahatMMA).^2);
        RiskJMA_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahatJMA).^2);
        RiskSAIC_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahatSAIC).^2);
        RiskSBIC_test(g,1) = sum((mu_test - X_test(:,1:M) * thetahatSBIC).^2);
        RiskLaMA_test(g,1) = sum((mu_test - X_test(:,1:M) *thetahatLaMA).^2);
    end
    %%  The relative training risk
    Risk = min(mean(Risk0));
    ARiskaic(r,1) = mean(Riskaic) / Risk;
    ARiskbic(r,1) = mean(Riskbic) / Risk;
    ARiskMMA(r,1) = mean(RiskMMA) / Risk;
    ARiskJMA(r,1) = mean(RiskJMA) / Risk;
    ARiskSAIC(r,1) = mean(RiskSAIC) / Risk;
    ARiskSBIC(r,1) = mean(RiskSBIC) / Risk;
    ARiskLaMA(r,1) = mean(RiskLaMA) / Risk;
    %%  The relative test risk
    Risk_test = min(mean(Risk0_test));
    ARiskaic_test(r,1) = mean(Riskaic_test) / Risk_test;
    ARiskbic_test(r,1) = mean(Riskbic_test) / Risk_test;
    ARiskMMA_test(r,1) = mean(RiskMMA_test) / Risk_test;
    ARiskJMA_test(r,1) = mean(RiskJMA_test) / Risk_test;
    ARiskSAIC_test(r,1) = mean(RiskSAIC_test) / Risk_test;
    ARiskSBIC_test(r,1) = mean(RiskSBIC_test) / Risk_test;
    ARiskLaMA_test(r,1) = mean(RiskLaMA_test) / Risk_test;

    fprintf('Progress: %d/%d\n', r, rl);
end
disp('--- In-sample loss ---');
result = round([ARiskaic, ARiskbic, ARiskSAIC, ARiskSBIC, ARiskJMA, ARiskMMA, ARiskLaMA], 4);
disp(['      AIC      BIC       SAIC      SBIC      JMA       MMA       LaMA']);
disp(result);

disp('--- Out-of-sample loss ---');
result_test = round([ARiskaic_test, ARiskbic_test, ARiskSAIC_test, ARiskSBIC_test, ARiskJMA_test, ARiskMMA_test, ARiskLaMA_test], 4);
disp(['      AIC      BIC       SAIC      SBIC      JMA       MMA       LaMA']);
disp(result_test);



coef = 0.1:0.1:0.9;
methods = {'AIC','BIC','SAIC','SBIC','JMA','MMA','LaMA'};
train_errors = result';
test_errors = result_test';
figure('Units','normalized','Position',[0.2 0.1 0.4 0.8])
t = tiledlayout(2, 1);
t.TileSpacing = 'tight'; 
t.Padding = 'compact';   
nexttile
h1 = plot(coef, train_errors(1,:),'--', 'LineWidth', 2, 'DisplayName', methods{1});
hold on
h2 = plot(coef, train_errors(2,:),':', 'LineWidth', 2, 'DisplayName', methods{2});
hold on
h3 = plot(coef, train_errors(3,:),'--', 'LineWidth', 2, 'DisplayName', methods{3});
hold on
h4 = plot(coef, train_errors(4,:),':', 'LineWidth', 2, 'DisplayName', methods{4});
hold on
h5 = plot(coef, train_errors(5,:),'-.', 'LineWidth', 2, 'DisplayName', methods{5});
hold on
h6 = plot(coef, train_errors(6,:),'r--', 'LineWidth', 2, 'DisplayName', methods{6});
hold on
h9 = plot(coef, train_errors(7,:),'k--square', 'LineWidth', 2, 'DisplayName', methods{7});
title('Relative in-sample loss',  'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
ylabel('Loss', 'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
xlabel('R^2', 'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
ylim([0.8 1.8])
xlim([0.1 0.9])
grid on
set(gca, 'LineWidth', 1, 'FontSize', 12, 'FontName', 'Times New Roman', 'Box', 'off')
lgd2 = legend('Location', 'southoutside','FontSize', 10);
legend('boxoff')
lgd2.NumColumns = 4;
nexttile
plot(coef, test_errors(1,:),'--', 'LineWidth', 2, 'DisplayName', methods{1})
hold on
plot(coef, test_errors(2,:),':', 'LineWidth', 2, 'DisplayName', methods{2})
hold on
plot(coef, test_errors(3,:),'--', 'LineWidth', 2, 'DisplayName', methods{3})
hold on
plot(coef, test_errors(4,:),':', 'LineWidth', 2, 'DisplayName', methods{4})
hold on
plot(coef, test_errors(5,:),'-.', 'LineWidth', 2, 'DisplayName', methods{5})
hold on
plot(coef, test_errors(6,:),'r--', 'LineWidth', 2, 'DisplayName', methods{6})
hold on
plot(coef, test_errors(7,:),'k--square', 'LineWidth', 2, 'DisplayName', methods{7})
title('Relative  out-of-sample loss', 'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
ylabel('Loss', 'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
xlabel('R^2', 'FontSize', 10,'LineWidth',3, 'FontName', 'Times New Roman')
ylim([0.8 1.8])
xlim([0.1 0.9])
grid on
set(gca, 'LineWidth', 1, 'FontSize', 12, 'FontName', 'Times New Roman', 'Box', 'off')
sgtitle(['n=',num2str(n),'  \alpha=',num2str(alpha),'  M=',num2str(M)], 'FontSize', 17, 'FontName', 'Times New Roman')



disp(['n',num2str(n),'a',num2str(alpha),'M',num2str(M)])



