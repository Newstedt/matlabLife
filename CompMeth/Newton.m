function root = newton(fun,x0,tol)
%NEWTON numerically solves equation system by using Newton?s mehod
%
% ROOT = NEWTON(FUN,X0,TOL) find an approximate solution ROOT to
% the equation system FUN(X)=0 f by using Newton?s method with XO
% as starting point and TOL as the stopping tolerance (that is,
% NORM(FUN(ROOT)) <= TOL). FUN is a function handle. The function
% [Y,J]=FUN(X) takes a vector argument X and returns the function
% value Y and the Jacobian J.
maxiter=100;
root=x0;
[y, J]=fun(root);
iter=0;

while norm(y)>tol && iter<maxiter
    
    root=root - J\y;
    [y, J]=fun(root);
    iter=iter + 1;
  
end
    
    