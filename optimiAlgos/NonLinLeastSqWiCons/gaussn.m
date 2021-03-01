function [x,code,n,X,alphas] = gaussn(f,X0,tol,maxIter,c1,aMin,params)
%GAUSSN Solves optimization problem using Gauss-Newton algorithm Armijo
%       linesearch.
%
%   IN:
%   f       - problem specific funciton
%   X0      - initial guess
%   tol     - minimum residual tolerance
%   maxIter - maximum iterations allowed before convergance
%   c1      - Armijo constant
%   aMin    - minimum alpha allowed before terminating
%   params  - cell array with any extra parameters required by f
%
%   OUT:
%   x       - final approximation of the solution
%   code    - code indicating success or failure
%   n       - number of consumed iterations
%   X       - succesive iterates of x
%   alphas  - a vector with the used step lengths
%
%   Minimum working example:
%   
%   %Define in-parameters
%   maxIter=50; tol=1e-8; c1=0.1; aMin=1e-3;
%
%   p = ellipse_data(Image(1),Image(2)); %Get datapoints
%   fun = @(y) Residual_Function(y,p); %Define function
%
%   >> [x,code,n,X,alphas] = gaussn(f,X0,tol,maxIter,c1,aMin,{})

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
%   2018-12-05: First version.

n = 0; % setup n
X(:,n+1) = X0; % assign starting guess for x
r = 0; % set r to ensure start of while loop
Fp0 = inf; % set Fp0 to ensure start of while loop

% loop while convergence criteria is not fulfilled as long as maximum 
% iterations has not been exceeded.
while (n < maxIter && norm(Fp0) >= tol*(1+norm(r)))
    n=n+1; % step
    
    % call function f to evaluate residuals and Jacobian for current x
    [r,J] = feval(f,X(:,n),params{:}); 
    
    Pk = (J'*J)\J'*(-r); % compute decent direction
    Fp0 = r'*J*Pk; % compute gradient
    
    % define objective function for calling linesearch
    f2 = @(x) 0.5*norm(feval(f,x,params{:}))^2;
    
    % call linesearch for current x
    alphas(n)= linesearch(f2,X(:,end),Pk,aMin,1/2*r'*r,Fp0,c1);
    X(:,n+1) = X(:,n) + alphas(n)*Pk; %save current x to X
end

x = X(:,end); % save final x for returning

% assign code reveiling success or failure/type of failure
if (alphas(end)<=aMin)
    code = -2;
elseif(n < maxIter)
    code = 0;
else
    code = -1;
end
