function [ranks, i] = pagerank(HT, alpha, eps)
%PAGERANK - Accepts a transposed sparse hyperlink matrix HT and performes 
%           PageRank on it. Also takes alpha (in Google matrix) and eps
%           (termination limit).          
%
%   MINIMAL WORKING EXAMPLE: We have a sparse hyperlink matrix AT and want
%                            to perform PageRank on it with alpha = 0.8 and
%                            eps = 1e-7.
%
%   >> [ranks, i] = pagerank(AT,0.8,1e-7)
%
%   ==> Score matrix: ranks
%       Number of iterations: i

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-12: Initial version .
%
% Function code starts here...

%pre-define initial pi with randomized elements
n = length(HT); %define the length of HT
pi = rand(n,1); %randomly generate a pi
piSum = sum(pi); %compute pi's sum
pi = pi ./ piSum; %all pi's elements by its sum to normalize it->sum(pi)=1
norm_pi = inf; %set the termination variable to inf 

%define a-matrix, holding ones for all-zero columns of HT and zeros else
a = zeros(n,1);
a((sum(HT,1) == 0)) = 1;

%define e-vector of all ones
e = ones(n,1);

%make all ingoing vectors sparse to preserve sparseness of HT
a = sparse(a);
e = sparse(e);
pi = sparse(pi);

% define/set iteration-counter
i = 0;

while norm_pi > eps %check terminating condition
    
    %save pi^(k)
    pi_prev = pi;
    
    %calculate pi^(k+1)
    pi = alpha*HT*pi...
         + (alpha/n)*e*(a'*pi)...
         + (1-alpha)/n*e*(e'*pi);
    
    %calculate beta, made for faster convergence
    beta = norm(pi_prev,1) - norm(pi,1);
    
    %apply beta on pi^(k+1)
    pi = pi + beta*(1/n);
    
    %calculate epsilon
    norm_pi = norm(pi - pi_prev);
    
    %count iteration
    i = i + 1;
end

%set ranks equal to pi for returning final result
ranks = pi;
