function root = Newton(fun,x0,tol)
%NEWTON solves equation system by using Newton?s mehod
%
% ROOT = NEWTON(FUN,X0,TOL) finds an approximate solution
% ROOT to the equation system FUN(X)=0 by using Newton?s
% method with XO as starting point and TOL as the stopping
% tolerance (that is, ABS(FUN(ROOT)) < TOL). FUN is a function
% handle. The function [Y,J]=FUN(X) takes a vector argument X
% and returns the function value Y and the Jacobian J. 

root=x0(:);
[y, j]=fun(root);
iter = 0;
maxiter=100;

while norm(y)>tol
    
    s = -j\y;
    root=root+s;
    iter = iter + 1;
    [y j]=fun(root);
    
end

