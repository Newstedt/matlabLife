function [x,code,n,X,alphas]=gaussn_niclas_damped(prob,x0,tol,maxIter,mu,alphaMin,params)
%GAUSSN Damped Gauss-Newton method, Goldstein-Armijo linesearch.
%
%[x,code,n,X,alphas]=gaussn(prob,x0,tol,maxIter,mu,alphaMin,params)
%prob     - name of problem.
%x0       - starting approximation.
%tol      - convergence tolerance.
%maxiter  - maximum number of iterations.
%alphaMin - shortest accepted step length.
%params   - cell array with additional arguments to residual and jacobian funs.
%x        - optimal solution.
%code     - error code
%              0 - OK
%             -1 - Too many iterations
%             -2 - No acceptable reduction found in line search.
%n        - number of consumed iterations.
%X        - iteration trace. X(:,i+1) is the ith iteration.
%alphas   - vector of accepted step lengths. alphas(i) is the ith step length.

% v1.0  1999-11-09. Niclas Borlin, niclas@cs.umu.se.
% v1.1  2002-05-07. Updated for the '02 course in Geometrical Image Analysis.
% v1.2  2003-05-19. Updated with explicit weights the '03 course in
%                   Geometrical Image Analysis.
% v1.3  2018-11-23. Added mu as a parameter.
