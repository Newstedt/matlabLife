function [insideInt] = intFun(u,Sigma,nu,s)

for j = 1:length(s)
    [z,w] = zFun(u,nu,s(j));
    phi = phiFun(Sigma,nu,z);
    
    wProd = 1;
    for i = 1:length(nu)
        wProd = wProd*(1/w(i));
    end

    insideInt(j) = phi*wProd;
end
