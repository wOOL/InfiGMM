function [X,Z,Mu,Sigma,Pi] = genData(N, D, K, Alpha)

    Mu = rand(K,D) * 30;

    Sigma = zeros(D,D,K);

    for k = 1:K
        Sigma(:,:,k) = wishrnd(eye(D), D+1);
    end

    %Pi = sampleDirichlet(ones(K,1),1)';
    
    Pi = stickBreaking(Alpha, K);

    [~, Z] = histc(rand(1,N),[0;cumsum(Pi)]);

    X = zeros(N,D);

    for k=1:K
        c = (Z == k);
        n = sum(c);
        [~,s,v] = svd(Sigma(:,:,k));
        X(c,:) = randn(n,D) * sqrt(s)*v' + repmat(Mu(k,:), n, 1);
    end
    
end

function Pi = stickBreaking(Alpha, K)

    V = betarnd(1,Alpha,[K,1]);
    V(end) = 1;
    rV = cumprod(1-V);
    Pi = [V(1);V(2:end).*rV(1:end-1)];
    
end 

% function Pi = sampleDirichlet(alpha, N)
% 
%     K = length(alpha);
%     Pi = zeros(N, K);
%     scale = 1;
%     for k=1:K
%       Pi(:,k) = gamrnd(alpha(k), scale, N, 1);
%     end
%     S = sum(Pi,2); 
%     Pi = Pi ./ repmat(S, 1, K);
% end