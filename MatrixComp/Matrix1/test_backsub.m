% TEST_BACKSUBST - Script for testing backsubst.m on randomly generated
%                  upper triangular matrix (of random size).
%
%   MINIMAL WORKING EXAMPLE: 
%   ">> test_backsubst" to see results in relevant test.

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

clear all
N = 50; %Define that random integers should be taken between 0 and N
n = 5; %Define size of nxn matrix

%Construct the randomised diagonal of the matrix random nxn upper
%triangular matrix
diagRand = randi(N,n);
diagRow = diagRand(1,:);

%Construct the randomised nxn upper triangular matrix with random integers
%between 0 and N
A = triu(randi(N,n),1) + diag(diagRow,0);
b = 10*rand(n,1); %Construct a randomised vector of matching size (nx1)

x = backsubst(A, b); %Solve system using backsubst
x2 = A\b; %Solve the system using backslash operator

normX = norm(x-x2) %Computes the norm of the difference between the methods

