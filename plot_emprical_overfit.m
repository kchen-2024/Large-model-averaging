clear all
n = 100;
p = 1000;
K = 50;
theta = ((1:p).^(-0.5));
theta =theta';
MM=[1:10:200,200:20:500];
for r = 1:length(MM)
    M=MM(r);
    s = 1:1:M;
    for g = 1:K
        %% Generating training data
        Xt = mvnrnd(zeros(1,p), eye(p), n); 
        Xt(:,1) = 1;  
        mu = Xt * theta;
        sigma=1;  
        e = randn(n,1)*sqrt(sigma);
        y = mu + e;
        %% Generating test data
        N_test=n;
        X_test = randn(N_test, p-1);
        X_test=[ones(N_test,1),X_test];
        sig2_test=X_test(:,2).^2;
        e_test=sig2_test.*randn(N_test, 1);
        mu_test=X_test*theta;
        Y_test=mu_test+e_test;
        %% Least squares estimation
        thetahat = zeros(p, M);
        for k = 1:M
            Xk = Xt(:,1:s(k));
            thetahat(1:s(k), k) = pinv (Xk'*Xk) * Xk'*y; 
        end
        w = ones(M,1)./M;     % equal weight
        Ln(r,g) = sum((mu - Xt*thetahat * w).^2)/n;
        test_Risk(r,g) = sum((mu_test - X_test*thetahat * w).^2)/ N_test;
        %% Delete the candidate model where kq=n
        if M>=n
            www = [ones(n-1,1);0;ones(M-n,1)]./(M-1);
            LLLn(r,g) = sum((mu - Xt*thetahat * www).^2)/n;
            tttest_Risk(r,g) = sum((mu_test - X_test*thetahat * www).^2)/ N_test;
        end
    end
    ELn(r)=mean(Ln(r,:));
    Etest_Risk(r)=mean(test_Risk(r,:));
    if M>=n
        ELLLn(r)=mean(LLLn(r,:));
        Etttest_Risk(r)=mean(tttest_Risk(r,:));
    else
        ELLLn(r)=ELn(r);
        Etttest_Risk(r)=Etest_Risk(r);
    end
    fprintf('Completion progress: %d/%d\n', r, length(MM));
end
%% Remove extreme values greater than 50
for r = 1:length(MM)
  tR=test_Risk(r,:);
        tR(tR>50)=[];
         EtR(r)=mean(tR);
end

%% Figure: Remove extreme values greater than 50
figure
plot(MM,ELn,'LineWidth',3)
hold on
plot(MM,EtR,'LineWidth',3)
xline(100,'--','LineWidth',2)
xlabel('M','LineWidth',5, 'FontName', 'Times New Roman');
ylabel('square error loss','LineWidth',5, 'FontName', 'Times New Roman');
lgd=legend('training squared error loss $L_n(\omega)$','test squared error loss $L_0(\omega)$','M=n', 'Interpreter', 'latex', 'FontName', 'Times New Roman' );
lgd.FontSize = 25;
lgd.FontName = 'Times New Roman';
set(lgd,'Box','off')
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 30);
set(gca,'Box','off')
ylim([0,12])

%% Figure: Delete the candidate model where kq=n
figure
plot(MM,ELLLn,'LineWidth',3)
hold on
plot(MM, Etttest_Risk,'LineWidth',3)
xline(100,'--','LineWidth',2)
xlabel('M','LineWidth',5, 'FontName', 'Times New Roman');
ylabel('square error loss','LineWidth',5, 'FontName', 'Times New Roman');
lgd=legend('training squared error loss $L_n(\omega)$','test squared error loss $L_0(\omega)$','M=n', 'Interpreter', 'latex', 'FontName', 'Times New Roman' );
lgd.FontSize = 25;
lgd.FontName = 'Times New Roman';
set(lgd,'Box','off')
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 30);
set(gca,'Box','off')
ylim([0,12])






