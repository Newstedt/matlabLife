%TESTMYPCG - Script that solves the questions 4.4.2-4.4.4 in Assignment 5.
%
%   MINIMAL WORKING EXAMPLE: 
%
%   >> testMyPCG;

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...

%% Performs tests in question 4.4.2
N = 100:100:600; %define N
lambda = [1,10,100]; %define lambda
tau = [1e-1 1e-2 1e-3 1e-4]; %define drop tolerance
tol = 1e-6; %define tolerance

%loop over all tau and lambda => 12 tables
for i = 1:length(tau)
    for j = 1:length(lambda)
        %print header
        fprintf('tau = %4.4f, lambda = %d \n',tau(i),lambda(j))
        
        %call myPCG for relevant lambda and tau
        myPCG(N,lambda(j),tau(i),tol);
        fprintf('\n \n') %print new line
    end
end

%% Time plots, i.e. tests in question 4.4.3
aux = cell(length(tau),1); %pre-define zero array for aux

%Loop over all tau to perform pcg using myPCG
for i = 1:length(tau)
   fprintf('tau = %f, lambda = %d \n',tau(i),lambda(j))
   aux{i} = myPCG(N,lambda(3),tau(i),tol); 
   fprintf('\n \n')
end

%plot the factorization time as a function of the drop tolerance
for i = 1:length(N)
    figure
    hold on
    for j = 1:length(tau)
        plotVar1(j) = aux{j}(i,5);
        plotVar2(j) = aux{j}(i,7);
        plotVar3(j) = aux{j}(i,8);
    end
    plot(tau,plotVar1,'linewidth',1.4)
    plot(tau,plotVar2,'linewidth',1.4)
    plot(tau,plotVar3,'linewidth',1.4)
    title(sprintf('N = %d', N(i)))
    xlabel('\tau')
    ylabel('time')
    legend('cholesky','pcg','tot')
    set(gca,'xscale','log')
    
    %!!!!! NOTE TO SELF: tau = 1e-3 gives lowest total time !!!!!!!!
end
    


