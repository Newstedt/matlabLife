function [x, flag, relres, it, resvec] = Jacobi(A, b, tol, maxit, x0)

% Jacobi  Jacobi's method for solving linear systems
%
% CALL SEQUENCE:
%
%      [x, flag, relres, it, resvec] = Jacobi(A, b, tol, maxit, x0)
%
% INPUT:
%    A       an n by n matrix
%    b       a compatible right hand side
%    tol     return when the relative residual drops below tol
%    maxit   execute at most maxit iterations
%    x0      the initial guess
%
% OUTPUT:
%    x       the computed approximation
%    flag    return flag,
%                flag = 0, if the routine converged
%                flag = 1, if the routine failed to converge
%    it      the number of iterations completed
%    resvec  the residual history
%
% MINIMAL WORKING EXAMPLE: TestJacobi.m

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2017-10-10  Initial programming and testing

% Isolate the diagonal
D=diag(diag(A)); 

% Define N;
N=D-A;

% Compute initial residual
r0=b-A*x0;

% Initialize output
x=x0; resvec=zeros(maxit+1,1); resvec(1)=norm(r0);

% Compute new right hand side
f=D\b;

% Assume failure to converge
flag=1;

% Main loop
for it=1:maxit
    % Update x
    x=D\(N*x)+f;
    % Compute residual
    r=b-A*x;
    % Save the norm of the residual for analysis
    resvec(it+1)=norm(r);
    % Check for convergence
    if (resvec(it+1)<resvec(1)*tol)
        flag=0; break;
    end
end

% Shrink output array to miminal size
resvec=resvec(1:it+1);

% Compute final relative residual
relres=resvec(end)/resvec(1);
