function B = triinv(A)
% TRIINV - Function that takes an nxn upper triangular matrix A and returns
%          the inverse of A.
%
%   Input parameters:
%   A: nxn upper triangular matrix
%
%   Output parameters:
%   B: Inverse of A
%
%   MINIMAL WORKING EXAMPLE: Find inverse of upper triangular matrix
%   defined as A = [1 3 6; 0 5 2; 0 0 9];
%
%   B = triinv(A);

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

B = zeros(size(A)); %pre-allocate room for B by making a zero matrix (for speed)
I = eye(length(A)); %create the nxn identity matrix
n = length(A); %define number of rows in A to loop over

%perform backsubstitution for each column of B to ultimatelly find the 
%inverse of A
for i = 1:n
    B(:,i) = backsubst(A,I(:,i));
end
