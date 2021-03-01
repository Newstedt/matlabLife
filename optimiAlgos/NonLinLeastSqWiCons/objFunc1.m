% First try to implement objective function for log-maximum likelihod
% estimisation. Two dimensions.
%
% INPUT
% N    : Number of simulations
% nu   : Column vector consisting of dof parameters
% Sigma: Correlation matrix in column vector-form, rows stacked on each
%        other.
%
% OUTPUT
% MLfunc : Log-likelihood function
%
% Example: N = 1000; nu = [5 20]'; Sigma = [1 0.8 0.8 1];
% 
% >> MLfunc = objFunc(N,nu,Sigma); 
%
% Code starts here...

function [MLfunc] = objFunc1(u,dim,theta)

nu = theta(1:dim);
A = reshape(theta(dim+1:end),[],sqrt(length(theta(dim+1:end))))';
Sigma = A*A';
% Maximum log-Likelihood Function, Part 1
for i = 1:length(u)
    insideInt = @(s) intFun(u(i,:)',Sigma,nu,s);
    logIntVal(i) = log(integral(insideInt,0,1));
end

MLp1 = sum(logIntVal);

% Maximum log-Likelihood Function, Part 2
dubSum = 0;
for i = 1:length(u)
    for j = 1:length(nu)
        dubSum = dubSum + 0.5*(nu(j) + 1)*log(1 + tinv(u(i,j),nu(j))^2/nu(j));
    end
end
MLp2 = dubSum;

% Maximum log-Likelihood Function, Part 3
Sum = 0;
for j = 1:length(nu)
    Sum = Sum + ( 0.5*log(nu(j)*pi) + log(gamma(0.5*nu(j))/gamma(0.5*(nu(j)+1))) );
end

MLp3 = length(u)*Sum;

% Put Maximum log-Likelihood Function parts together
MLfunc = MLp1 + MLp2 + MLp3;
