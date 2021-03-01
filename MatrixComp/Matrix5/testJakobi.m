% TestJacobi MWE for Jacobi routine

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2017-10-10  Initial programming and testing

% Make the system
MakeSystem;

% Apply Jacobi's method
tic; [x, flag, relres, it, resvec]=Jacobi(A,b,1e-6,100,x0); toc;