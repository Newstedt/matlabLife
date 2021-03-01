function [x,code,n,X,dk,rhos,stepType] = levenberg(f,x0,tol,maxIter,d0,mu,eta,params)
%LEVENBERG Solves optimization problem using geometric 
%          Levenberg-Marquardt-Powell algorithm.
%
%   IN:
%   f       - problem specific funciton
%   X0      - initial guess
%   tol     - minimum convergence tolerance
%   maxIter - maximum iterations allowed before convergance
%   d0      - Initial trust region
%   mu, eta - the good/bad thresholds
%   params  - cell array with any extra parameters required by f
%
%   OUT:
%   x        - final approximation of the solution
%   code     - code indicating success or failure
%   n        - number of consumed iterations
%   X        - succesive iterates of x
%   dd       - a vector with trust regions
%   rhos     - a vector with gain ratios
%   stepType - vector indicating which type of steps that has been taken
%
%   Minimum working example:
%   
%   %Define in-parameters
%   maxIter=50; tol=1e-8; c1=0.1; aMin=1e-3; mu=0.25; eta=0.75; dk=0; 
%   rhos=0;
%
%   p = ellipse_data(Image(1),Image(2)); %Get datapoints
%   fun = @(y) Residual_Function(y,p); %Define function
%
%  [x,code,n,X,dk,rhos,stepType] = levenberg(f,x0,tol,maxIter,d0,mu,eta,{})

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
%   2018-12-05: First version.

n = 1; % set counter
dk = d0; % set initial trust region
X(:,n) = x0; % set initial guess

% call function f to evaluate residuals and Jacobian for initial x
[r,J] = feval(f,X,params{:}); 
Accepted(1) = 1; %pre-set accept-variable to true

% loop while convergence criteria is not fulfilled as long as maximum 
% iterations has not been exceeded.
while(norm(J*((J'*J)\(J'*(-r)))) >= tol*(1+norm(r)) && n < maxIter)
    
    %if last trust region was accepted, do...
    if(Accepted(end))
        % call function f to evaluate residuals and Jacobian for initial x
        [r,J] = feval(f,X(:,end),params{:});
        
        Pgn = (J'*J)\(J'*(-r)); % compute Gauss-Newton decent direction
        g = J'*r; % define g for computing cauchy point
        Pcp = -((g'*(J')*J*g)\(g'*g))*(g); % compute cauchy point
    end
    
    %if Gauss-Newton decent direction is less or equal to last trust region
    if (norm(Pgn)<=dk(end))
        P = Pgn; % set P to Gauss-Newton decent direction
        stepType(n) = 0; % set step type

    elseif(norm(Pcp) >= dk(end))
        P = dk(end)*Pcp/norm(Pcp);
        stepType(n) = 2; % set step type

    else
        % solve for intersection between Cauchy point and Pgn
        a = (Pcp - Pgn);
        b = Pgn;
        A = a'*a;
        B = a'*b+b'*a;
        C = b'*b;
        beta = (-sqrt(-4*A*C+4*A*dk(end)^2+B^2)-B)/(2*A);
        
        P = Pcp*beta+Pgn*(1-beta); %set P
        stepType(n) = 1; % set step type
    end
    
    tk = X(:,end)+ P; % calculate trial point
    
    % calculate gain ratio
    rtk = feval(f,tk,params{:}); 
    rhos(n) = (1/2*r'*r - 1/2*rtk'*rtk)/(-P'*J'*r-1/2*P'*J'*J*P);

    % check if we should accept current trust region and update it in 
    % accordance for next iteration
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
    n = n+1; % step
end
x = X(:,end); % save final x for returning 
if (maxIter == n)
    code = -1; % assign code
else
    code = 0; % assign code
    n=n-1;
end
