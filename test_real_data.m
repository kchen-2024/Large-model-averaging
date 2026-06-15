clear all
%% The data after sorting the importance of variables
data = readtable('ordered_mtcars.csv');
% data = readtable('ordered_UScrime.csv');
DATA = table2array(data(:, 1:end-1));    % Predictor variable
yy = data.y;   % Response variable

%% Size of training set
nn=[13;16;19;22;25;28];  % mtcars
% nn=[18;21;24;27;30;33];  % UScrime
G=10000;
for num=1:length(nn)
    Riskaic = zeros(G, 1);
    Riskbic = zeros(G, 1);
    RiskSAIC = zeros(G, 1);
    RiskSBIC = zeros(G, 1);
    RiskMMA = zeros(G, 1);
    RiskJMA = zeros(G, 1);
    RiskRMT = zeros(G, 1);
    for g=1:G
        %% Divide the training set and the test set
        M=size(DATA,2)+1;
        s=1:M;
        N=size(DATA,1);
        n=nn(num);
        A=randperm(N);
        B=A(1:n);
        C=A(n+1:N);
        Xtrain=DATA(B,:);
        ytrain=yy(B);
        Xtest=DATA(C,:);
        ytest=yy(C);
        %% Standardized training set and test set
        for i = 1 : size(DATA,2)
            mea = mean(Xtrain(:,i));
            st = std(Xtrain(:,i));
            if st > 1e-12
                Xtrain(:,i) = (Xtrain(:,i)-mea)/st;
                Xtest(:,i) = (Xtest(:,i)-mea)/st;
            else
                Xtrain(:,i) = 0;
                Xtest(:,i) = 0;
            end
        end
        ymea = mean(ytrain);
        yst = std(ytrain);
        ytrain = (ytrain-ymea)/yst;
        ytest = (ytest-ymea)/yst;

        Xtrain = [ ones(n,1) , Xtrain];
        Xtest = [ ones(N-n,1) , Xtest];
        Xt = Xtrain;
        X = Xt(:,1:M);
        y = ytrain;
        %% JMA
        ehat_jack = zeros(n,M);
        mu_tilde=zeros(n,M);
        for k = 1:M
            X_sub = X(:,1:k);
            Pmm= X_sub*inv(X_sub'*X_sub+eye(k)*1e-6)* X_sub';
            h_ii=diag(Pmm);
            Dmm=diag(1./max(1-h_ii,eps));  % Add a small tolerance to prevent division by zero
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
        % mm=round(3 * (n^(1/3)));
        % mm=M;
        % H = Xt(:,1:mm) *pinv (Xt(:,1:mm)'*Xt(:,1:mm)) * Xt(:,1:mm)';
        % e_large =(eye(n) - H) * y;
        % sigma = e_large' *e_large / (n - mm);

        thetahat = zeros(M,M);
        for k = 1:M
            Xk = X(:,1:k);
            thetahat(1:k, k) = inv(Xk' * Xk++eye(k)*1e-6) * (Xk' * y);
        end
        muhat = X * thetahat;
        SSR = sum((y - muhat).^2, 1);
        Y = repmat(y, 1, M);
        ehat = Y - muhat;
        a1 = ehat' * ehat;
        if rank(a1) < size(ehat,2)
            a1 = a1 + eye(M)*1e-10;
        end
        %% noise variance estimated by BIC-selected model
        bic = n * log(SSR)- n * log(n)+ s * log(n);
        [~, bicch] = min(bic);
        sigma = SSR(bicch) / (n - s(bicch));
        a2 = sigma * s(:);
        w = quadprog(a1, a2, [], [], Aeq, beq, lb, ub, w0, options);
        w = max(w, 0);
        w = w ./ sum(w);
        %% LaMA
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
        w_RMT = quadprog(2*(a1+n*Q_V+n*An), zeros(M, 1), [], [], Aeq, beq, lb, ub, w0, options);
        w_RMT = max(w_RMT, 0);
        w_RMT = w_RMT ./ sum(w_RMT);
        thetahatRMT = thetahat * w_RMT;
     
        %% AIC/BIC
        aic = n * log(SSR) + s * 2 - n * log(n);
        [~, aicch] = min(aic);
        wa = exp(-0.5 * aic);
        wa = wa / sum(wa);
        wb = exp(-0.5 * bic);
        wb = wb / sum(wb);
        %% Calculate the estimators of each method
        thetahataic = thetahat(:, aicch);
        thetahatbic = thetahat(:, bicch);
        thetahatMMA = thetahat * w;
        thetahatJMA = thetahat * w_jack;
        thetahatSAIC = thetahat * wa';
        thetahatSBIC = thetahat * wb';
        %% test error
        RiskMMA(g) = mean((ytest - Xtest*thetahatMMA).^2);
        RiskRMT(g) = mean((ytest - Xtest*thetahatRMT).^2);
        Riskaic(g) = mean((ytest  - Xtest * thetahataic).^2);
        Riskbic(g) = mean((ytest  - Xtest * thetahatbic).^2);
        RiskJMA(g) = mean((ytest  - Xtest * thetahatJMA).^2);
        RiskSAIC(g) = mean((ytest  - Xtest * thetahatSAIC).^2);
        RiskSBIC(g) = mean((ytest - Xtest * thetahatSBIC).^2);
    end
    me(num,:)=[mean(Riskaic),mean(Riskbic),mean(RiskSAIC),mean(RiskSBIC),mean(RiskMMA),mean(RiskJMA),mean(RiskRMT)];
    va(num,:)=[var(Riskaic),var(Riskbic),var(RiskSAIC),var(RiskSBIC),var(RiskMMA),var(RiskJMA),var(RiskRMT)];
end
%% Compare the mean and variance of the test errors
disp('   n          AIC       BIC       SAIC      SBIC      MMA       JMA       LaMA')
disp([nn,me])
disp([nn,va])

results=[nn,me;nn,va];
