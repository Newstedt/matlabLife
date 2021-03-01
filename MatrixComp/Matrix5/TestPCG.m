% TestPCG 

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2017-10-10  Initial programming and testing

% Make system
MakeSystem;

% Construct basic preconditoner
L=ichol(A);

% Set convergence tolerance and maximum number of iterations
tol=1e-6; maxit=100;

% Run the PCG algorithm
tic; [x, flag, relres, iter, resvec]=pcg(A,b,tol,maxit,L,L'); toc;

% Set options for regular, incomplete Cholesky with droptol 1e-4
opts=struct('type','ict','droptol',1e-4,'michol','off');

% Try the more advanced preconditioner
L=ichol(A,opts);

% Run the PCG algorithm again
tic; [x, flag, relres, iter, resvec]=pcg(A,b,tol,maxit,L,L'); toc;

