function [x,code,n,X,alphas]=gaussn_niclas_undamped(prob,x0,tol,maxIter,params)
%GAUSSN Undamped Gauss-Newton method.
%
%[x,code,n,X,alphas]=gaussn(prob,x0,tol,maxIter,params)
%prob     - name of problem.
%x0       - starting approximation.
%tol      - convergence tolerance.
%maxiter  - maximum number of iterations.
%params   - cell array with additional arguments to residual and jacobian funs.
%x        - optimal solution.
%code     - error code
%              0 - OK
%             -1 - Too many iterations
%             -2 - No acceptable reduction found in line search.
%n        - number of consumed iterations.
%X        - iteration trace. X(:,i+1) is the ith iteration.
%alphas   - vector of accepted step lengths. alphas(i) is the ith step length.

