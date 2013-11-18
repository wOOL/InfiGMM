function [gamma,mu_0,lambda,W,nu] = Mstep(X,Phi,alpha)

    D = size(X,2);
    T = size(Phi,2);
    mu_hat = mean(X);
    W_0 = .1 * D * cov(X);

    Ns = sum(Phi,1) + 1e-10;
    sigs = zeros(D,D,T);
    mus = Phi' * X ./ repmat(Ns',1,D);
    for t = 1:T
        diff0 = X - repmat(mus(t,:),size(X,1),1);
        diff1 = repmat(sqrt(Phi(:,t)),1,D) .* diff0;
        sigs(:,:,t) = diff1' * diff1;
    end

    gamma(:,1) = 1 + sum(Phi);
    gamma(:,2) = alpha + flipdim(cumsum(flipdim(sum(Phi),2)),2) - sum(Phi);
    lambda = Ns + 1;
    nu = Ns + D;
    mu_0 = (repmat(Ns',1,D) .* mus + repmat(mu_hat,size(Ns,2),1)) ./ repmat(Ns' + 1,1,D);
    W = zeros(D,D,T);
    for t = 1:T
        diff = mus(t,:) - mu_hat;
        W(:,:,t) = sigs(:,:,t) + Ns(t) * diff * diff' / (Ns(t)+1) + W_0;
    end
    
end