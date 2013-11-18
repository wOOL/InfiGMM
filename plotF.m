function plotF(data, z)

    S = ['d','o','x','+','.','v','^'];
    C = ['y','m','c','r','g','b','k'];
    for i = min(z):max(z)
        hold on;
        scatter(data(z==i,1),data(z==i,2),S(mod(i,length(S))+1),C(mod(i,length(C))+1));
    end
    legend;
    
end