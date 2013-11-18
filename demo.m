close all;
clear;
% Set generative parameters
rng(1943);
N = 300;
D = 2;
T = 3;
Alpha = 1;
% Generate data points
[X,Z,Mu,Sigma,Pi] = genData(N,D,T,Alpha);
% Plot
subplot(1,2,1);
title('True');
plotF(X,Z);
% Initialise
T = 50;
Phi = 1/T*ones(N,T);
alpha = rand*10;
% Run EM!
for i = 1:100
    [gamma,mu_0,lambda,W,nu] = Mstep(X,Phi,alpha);
    [Phi,alpha] = Estep(X,gamma,mu_0,lambda,W,nu);
end
% Plot result
subplot(1,2,2);
title('Estimate');
[~, Z_hat] = max(Phi,[],2);
Z = zeros(1,N);
UZ = sort(unique(Z_hat));
for i = 1:length(UZ)
    Z(Z_hat==UZ(i)) = i;
end
plotF(X,Z);