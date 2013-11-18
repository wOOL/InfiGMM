function [Phi,alpha] = Estep(X,gamma,mu_0,lambda,W,nu)

    N = size(X,1);
    T = size(nu,2);
    Eq_log_V = zeros(T,1);
    Eq_log_1_V = zeros(T,1);
    Eq_log_Z = zeros(T,1);
    S = zeros(N,T);
    log_X = zeros(N,T);
    for t = 1:T
        Eq_log_V(t) = psi(gamma(t,1)) - psi(gamma(t,1)+gamma(t,2));
        Eq_log_1_V(t) = psi(gamma(t,2)) - psi(gamma(t,1)+gamma(t,2));
        Eq_log_Z(t) = Eq_log_V(t) + sum(Eq_log_1_V(1:t-1));
        log_X(:,t) = normalWishart(X,mu_0(t,:),lambda(t),W(:,:,t),nu(t));
        S(:,t) = Eq_log_Z(t) + log_X(:,t);
    end
    Phi = exp(S);
    Phi = Phi ./ repmat(sum(Phi,2),1,T);
    alpha = T / (1 - sum(Eq_log_1_V(1:(end-1))));
    
end

function logprob = normalWishart(X,mu_0,lambda,W,nu)

    [N,D] = size(X);
    d = X - squeeze(repmat(mu_0,N,1));
    logprob = -0.5 * D * log(2*pi) - 0.5 * logdet(0.5*W) + ...
               0.5 * sum(psi(0.5 * (nu + 1 - (1:D)))) - 0.5 * D / lambda- ...
               sum((0.5 * nu * d / W).*d,2);
             % sum((0.5 * nu * d * inv(W)).*d,2);
    
end

function y = logdet(A)
    
    U = chol(A);
    y = 2*sum(log(diag(U)));
    
end