function [x,code,n,X,alphas] = gaussn(f,X0,tol,maxIter,c1,aMin,params)
n = 0;
X(:,n+1) = X0;
r = 0;
Fp0 = inf;

while (n < maxIter && norm(Fp0) >= tol*(1+norm(r)))
    n=n+1;
    [r,J] = feval(f,X(:,n),params{:});
    %Pk = -J\r;
    Pk = (J'*J)\J'*(-r);
    Fp0 = r'*J*Pk; %-r;
    f2 = @(x) 0.5*norm(feval(f,x,params{:}))^2;
    alphas(n)= linesearch(f2,X(:,end),Pk,aMin,1/2*r'*r,Fp0,c1);
    X(:,n+1) = X(:,n) + alphas(n)*Pk;
end

x = X(:,end);
% 
if (alphas(end) <= aMin)
    code = -2;
elseif(n < maxIter)
    code = 0;
else
    code = -1;
end

n = n-1;