function plotF(X, Z)

    S = ['d','o','x','+','.','v','^'];
    C = ['y','m','c','r','g','b','k'];
    for i = min(Z):max(Z)
        hold on;
        scatter(X(Z==i,1),X(Z==i,2),S(mod(i,length(S))+1),C(mod(i,length(C))+1));
    end
    legend;
    
end