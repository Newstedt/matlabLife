%TESTCG - Script that determines the smallest value of N, that will make 
%          the non-preconditioned CG algorithm require more than one second
%          to solve a linear system of size A.
%
%   MINIMAL WORKING EXAMPLE: 
%
%   >> testCG;

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...
N = 1400:10:1600; %define range for N
M = N-1; %define M
lambda=1; %define chritical parameter
tol = 1e-6; %set tolerance
maxit = 100; %set maximum iterations before termination
resvec = cell(length(N),1); %pre-define empty cell array for residuals
time = zeros(length(N),1); %pre-define zero vector for storing time

%loop over all N
for i = 1:length(N)
    A=Heat(N(i),lambda); %construct heat equation matrix
    seed = 2019; %(omg future seed)
    rng(seed);
    b=rand(M(i)^2,1); %construct randomized b vector
    tic %start timing
    [x, flag, relres, iter, resvec]=pcg(A,b,tol,maxit); %call pcg for systm
    time(i) = toc; %collect time spent
end

N(find(time>1,1)) %print N-limit asked for. 1440 for example (one occasion)

