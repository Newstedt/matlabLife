% TEST_TRIINV - Script for testing triinv.m on four randomly generated
%               upper triangular matrix of size 3x3, 4x4, 5x5 and 6x6.
%
%   MINIMAL WORKING EXAMPLE: 
%   ">> test_triinv" to see results in relevant test.

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...
clear all
k = 1;
for i = 3:6
n = i; %Define size of nxn matrix

%Construct the randomised nxn upper triangular matrix with random numbers
%between 0 and 100
A = triu(100*rand(n),1) + diag(100*rand(n,1),0);

B = triinv(A); %Compute inverse of A by using triinv(A)

%Check correctness by computing the Frobenius matrix norm of AB - I
Frob(k) = norm(A*B-eye(n),'fro');
k = k + 1;
end

Frob = Frob'
