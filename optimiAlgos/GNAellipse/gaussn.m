function [x,code,n,X,alphas] = gaussn(f,x0,tol,maxIter,c1,aMin,params)

n = 0;
X(:,1) = x0;
code = 0;

%Set initial values for r and Fp0 to ensure we start the loop
Fp0 = inf; r = 0;

while n < maxIter && norm(Fp0) >= tol*(1+norm(r))
    
    n = n + 1;
    [r,J] = feval(f,X(:,n),params{:});
    pk = -J\r;
    Fp0 = r'*J*pk;
    F0 = 0.5*r'*r;
    f2 = @(x) 0.5*norm(feval(f,x,params{:}))^2;
    alphas(n) = linesearch(f2,X(:,n),pk,aMin,F0,Fp0,c1);
    if alphas(n) < aMin
        code = 2;
    end
    X(:,n+1) = X(:,n) + alphas(n)*pk;
    
end

if n >= maxIter && code ~= 2
    code = 1;
end

x = X(:,end);
    
    