function [X] = myRandTsampNuG(N,nuG,covMat)
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

X = zeros(N,length(covMat));

for j = 1:N
    S = rand;
    Z = mvnrnd(zeros(1,length(covMat)),covMat,1);
    W = sqrt(nuG/chi2inv(S,nuG));
    X(j,:) = W*Z;
end