function [x,code,n,X,dk,rhos,stepType] = levenberg(f,x0,tol,maxIter,d0,mu,eta,params)
n = 1;
dk = d0;
X(:,n) = x0;
[r,J] = feval(f,X,params{:});
Accepted(1) = 1;
rhos(1)=0;
while(norm(J*((J'*J)\(J'*(-r)))) >= tol*(1+norm(r)) && n < maxIter)
    if(Accepted(end))
        [r,J] = feval(f,X(:,end),params{:});
        Pgn = (J'*J)\(J'*(-r));
        g = J'*r;
        Pcp = -((g'*(J')*J*g)\(g'*g))*(g);
    end
    
    if (norm(Pgn)<=dk(end))
        P = Pgn;
    elseif(norm(Pcp) >= dk(end))
        P = dk(end)*Pcp/norm(Pcp);
    else
        a = (Pcp - Pgn);
        b = Pgn;
        A = a'*a;
        B = a'*b+b'*a;
        C = b'*b;
        beta = (-sqrt(-4*A*C+4*A*dk(end)^2+B^2)-B)/(2*A);
        P = Pcp*beta+Pgn*(1-beta);
    end
    tk = X(:,end)+ P;
    rtk = feval(f,tk,params{:});
    rhos(n) = (1/2*r'*r - 1/2*rtk'*rtk)/(-P'*J'*r-1/2*P'*J'*J*P);

    if (rhos(end) < mu)
        X(:,end+1)=X(:,n);
        dk(n+1) = dk(end)/2;
        Accepted(n+1) = 0;
    elseif(rhos(end) < eta)
        X(:,end+1) = tk;
        dk(n+1)=dk(end);
        Accepted(n+1) = 1;
    else
        X(:,end+1) = tk;
        dk(n+1) = 2*dk(end);
        Accepted(n+1) = 1;
    end
    n=n+1;
end
x = X(:,end);
if (maxIter == n)
    code = -1;
else
    code = 0;
    n=n-1;
end
stepType=0;
