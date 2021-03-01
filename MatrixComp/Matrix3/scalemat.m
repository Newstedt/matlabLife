function B = scalemat(A)
%SCALEMAT - Takes a matrix A, with binary elements (0 or 1) and re-scales
%           it such that each column sum to 1. Returns sparse matrix B.
%
%   MINIMAL WORKING EXAMPLE: Rescale matrix,
%
%           A= [1 0 1 1;
%               0 0 1 1; 
%               1 0 0 1;
%               0 0 1 0],
%           
%           such that each column of the rescaled matrix B sum to 1.
%
%   >> A = [1 0 1 1; 0 0 1 1; 1 0 0 1; 0 0 1 0]; %define A
%   >> B = scalemat(A); %re-scale A by using "scalemat"
%   >> full(B) %visualize your result by computing the full matrix

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-08: Initial version .
%
% Function code starts here...

%Make A sparse
A = sparse(A); 

%Define vector of scaling factors for each column. Put them in a diagonal-
%sparse matrix.
scale = diag(sparse(1./sum(A,1)));

%Scale A and set this as B.
B = A*scale;





