function [x,code,n,X,dd,rhos,stepType] = levenberg(f,x0,tol,maxIter,d0,mu,eta,params)

X(:,1) = x0;
k = 0;
n = 1;
dd(1) = d0;
[r,J] = feval(f,X,params{:});
accepted(1) = 1;
rhos(1) = 0;

%----------------------- 2 ---------------------------------------
while(norm(J*((J'*J)\(J'*(-r)))) >= tol*(1+norm(r)) && n < maxIter)
    %------------------------- 2.1 -------------------------------
    if accepted(end) == 1
        [r,J] = feval(f,X(:,n),params{:});
        pGn = (J'*J)\(J'*(-r));
        g = J'*r;
        pCp = -((g'*(J')*J*g)\(g'*g))*(g);
    end
    %------------------------- 2.2 -------------------------------
    if norm(pGn) <= dd(end)
        pk = pGn;
        stepType(n) = 0;
    elseif norm(pCp) >= dd(end)
        pk = dd(end)*pCp/norm(pCp);
        stepType(n) = 2;
    else
        a = (pCp - pGn);
        b = pGn;
        A = a'*a;
        B = a'*b+b'*a;
        C = b'*b;
        beta = (-sqrt(-4*A*C+4*A*dd(end)^2+B^2)-B)/(2*A);
        pk = pCp*beta+pGn*(1-beta);
        stepType(n) = 1;
    end
    %------------------------- 2.3 -------------------------------
    t = X(:,end) + pk;
    %------------------------- 2.4 -------------------------------
    rhos(n) = (feval(f,X(:,end))-feval(f,X(:,end)+pk))/(-pk'*J'*r-0.5*pk'*J'*J*pk);
    %------------------------- 2.5 -------------------------------
    if rhos(end) < mu
        X(:,n+1) = X(:,n);
        dd(n+1) = dd(n)/2;
        accepted(n+1) = 0;
    elseif rhos(n) < eta
        X(:,n+1) = t;
        dd(n+1) = dd(n);
        accepted(n+1) = 1;
    else
        X(:,n+1) = t;
        dd(n+1) = 2*dd(n);
        accepted(n+1) = 1;
    end
    n = n + 1;
end

if n > maxIter
   code = -1; 
else
   code = 0
   n = n-1;
end

x = X(:,end);


