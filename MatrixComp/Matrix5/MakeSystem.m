% MakeSystem  Generates a system corresponding to 2D heat equation

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2017-10-10  Initial programming and testing

% Dimension, M is the number of internal grid points in each direction
N=2000; M=N-1;

% Set the critical parameter \ttlambda
lambda=1;
 
% Matrix
A=Heat(N,lambda);

% Seed
seed=2017;

% Initialize generator
rng(seed);

% Random right hand side
b=rand(M^2,1);

% Zero Initial guess
x0=zeros(M^2,1);