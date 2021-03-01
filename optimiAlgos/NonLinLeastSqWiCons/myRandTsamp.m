function [X] = myRandTsamp(N,nu,Sigma)
%   Function that produces randomized values from multivariate
%   t-distribution, with different degrees of freedom in its univariate
%   distributions.
%
%   INPUT 
%   N    : Number of randomized values from distribution
%   nu   : Vector of different degrees of freedom for the correlated
%          distributions
%   Sigma: Correlation matrix
%
%   OUTPUT
%   X    : Vector of random values
%
%   Code starts here...

X = zeros(N,length(nu));

for j = 1:N
    S = rand;
    covMat = corr2cov([1 1],Sigma);
    Z = mvnrnd([0 0],covMat,1);
    for i = 1:length(nu)
        W = sqrt(nu(i)/chi2inv(S,nu(i)));
        X(j,i) = W*Z(i);
    end
end