function [aux] = myPCG(N,lambda,tau,tol)
%MYPCG - Constructs a linear system based on given parameters, and then
%        solves it using pcg. Prints result for all N given
%
%   MINIMAL WORKING EXAMPLE: Construct and solve systems with size
%                            N = 10:10:100, with critical parameter 
%                            lambda=1, drop tolerance tau = 1e-3 and
%                            tolerance tol = 1e-6:
%
%   >> aux = myPCG([10:10:100],1,1e-3,1e-6);

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...

M = N-1; %define M

aux = zeros(length(N),8); %pre-define zero array for aux
p = 0; %pre-define termination variable for onZeroL/nonZeroA condition
maxit = 100; %set maximum iterations
 
% Loop over all N
for i = 1:length(N)
    
    A=Heat(N(i),lambda); %Construct A using Heat.m
    nonZeroA = length(nonzeros(A)); %Find non-zero elements in A
    
    %set options for ichol
    opts=struct('type','ict','droptol',tau,'michol','off'); 
    tic %start timing for ichol
    L = ichol(A,opts); %construct L by using ichol
    time = toc; %stop timing, collect time spent
    nonZeroL = length(nonzeros(L)); %find non-zero elements in L
    
    seed = 2019; %(omg future seed)
    rng(seed); %set seed
    b=rand(M(i)^2,1); %get randomized b vector
    
    tic %start timing for pcg
    %solve system using pcg
    [x, flag, relres, iter, resvec]=pcg(A,b,tol,maxit,L,L');
    timeSolve = toc; %stop timing for pcg, collect time spent
    
    timeSum = time+timeSolve; %sum the total time spent
    
    %set values in aux for printing
    aux(i,:) = [N(i), M(i)^2, nonZeroA, nonZeroL,...
                time, flag, timeSolve, timeSum];
    
    %check if ratio has exceeded 10 and if it is the first time
    if nonZeroL/nonZeroA > 10 && p == 0
        exceed = N(i); %N for which ratio exceeds 10
        p = 1; %set first time check-variable 
    end
    
end

%define header string
text1 = ['       N     M^2   nnz(A)    nnz(L)     time  '...
         '\t flag  timeSolve   timeSum\n'];
fprintf(text1); %print header string

%below is code for latex insertion
% text2 = ['N & M^2 & nnz(A) & nnz(L) & time & flag & timeSolve'...
%          '& timeSum \\\\ \\hline \n'];
% fprintf(text2);

for i = 1:length(N)
    %set row string
    text3 = ['%8d %8d %8d %10d %10.2e   %d %10.2e %10.2e\n'];
    fprintf(text3,aux(i,:)) %print row string
    
%below is code for latex insertion
%     text4 = ['%8d&%10.2e&%10.2e&%10.2e&%10.2e&%d&%10.2e&%10.2e' ...
%         '\\\\ \\hline \n'];
%     fprintf(text4,aux(i,:))

end


