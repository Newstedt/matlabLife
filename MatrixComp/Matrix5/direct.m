%DIRECT - Script that runs tests on direct solvers and prints the result
%         from questions in 4.1 of the assignment
%
%   MINIMAL WORKING EXAMPLE: 
%
%   >> direct;

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...
N = 10:10:400; %define range of N
M = N-1; %from N, define M

lambda=1; %set the critical parameter lambda
time = zeros(length(N),1); %pre-define zero vector for time
nonZeroA = zeros(length(N),1); %pre-define zero vetor for nonZeroA
nonZeroL = zeros(length(N),1); %pre-define zero vetor for nonZeroL
aux = zeros(length(N),6); %pre-define zero vetor for aux
p = 0; %set criterion variable
 
%loop over all N
for i = 1:length(N)
    A=Heat(N(i),lambda); %construct heat equation matrix
    nonZeroA = length(nonzeros(A)); %store non-zero elements in A
    tic %start timing
    L = chol(A); %perform Cholesky factorization on A
    time = toc; %store time spent on chol(A)
    nonZeroL = length(nonzeros(L)); %store non-zero elements in L
    b = ones(length(L),1); %set b to vector of all ones
    tic %start timing
    x = L\b; %solve system with cholesky-L
    timeSolve = toc; %store time for solving system
    
    %set variables for printing
    aux(i,:) = [N(i), M(i)^2, nonZeroA, nonZeroL, time, timeSolve];
    
    %check if ratio is larger than 10 (for the first time (--> p))
    if nonZeroL/nonZeroA > 10 && p == 0
        exceed = N(i); %N for which ratio exceeds 10
        p = 1; %set criterion variable
    end
    
end

%print table
fprintf('       N      M^2   nnz(A)     nnz(L)     time     timeSolve \n');
for i = 1:length(N)
    fprintf('%8d %8d %8d %10d %10.2e %10.2e \n',aux(i,:))
end

timeExc = find(aux(:,5) > 1,1); %time for which factorization time exceed 1
NtimeExc = N(timeExc); %N for which factorization time exceed 1
