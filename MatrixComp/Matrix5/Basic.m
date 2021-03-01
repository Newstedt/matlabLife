%BASIC - Script that runs tests on Jacobi's method and prints the result
%        from questions in 4.2 of the assignment
%
%   MINIMAL WORKING EXAMPLE: 
%
%   >> Basic;

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...
N = 100:100:600; %define N
M = N-1; %define M
lambda=1; %set critical parameter lambda
tol = 1e-6; %set tolerance
maxit = 100; %set maximum iterations before definite termination
resvec = cell(length(N),1); %pre-define empty cell array for residuals

%loop over all N
for i = 1:length(N)
    A=Heat(N(i),lambda); %construct heat equation matrix A
    seed = 2019; %(omg future seed)
    rng(seed);
    b=rand(M(i)^2,1); %construct randomized b vector for system
    x0 = zeros(M(i)^2,1); %set initial condition
    
    %solve system using Jacobi's method
    [x, flag, relres, it, resvec{i}] = Jacobi(A, b, tol, maxit, x0);
end

%% Plot residual history
figure
hold on
for i = 1:length(N)
    plot(resvec{i})
    xlabel('iter')
    ylabel('resVec')
    legendInfo{i} = ['N = ' num2str(N(i))];
end
legend(legendInfo);

%% Linear convergence
figure
hold on 
for i = 1:length(N)
%     plot(resvec{i}(1:end-1)./resvec{i}(2:end))
    plot(resvec{i}(2:end)./resvec{i}(1:end-1))
    xlabel('iter')
    ylabel('err_{i+1}/err_i')
    legendInfo{i} = ['N = ' num2str(N(i))];
end
legend(legendInfo);

%% N independence of Jacobi method
figure
hold on
for i = 1:length(N)
    plot(log(resvec{i}/N(i)))
    xlabel('iter')
    ylabel('log(resVec/N)')
    legendInfo{i} = ['N = ' num2str(N(i))];
end
legend(legendInfo);

