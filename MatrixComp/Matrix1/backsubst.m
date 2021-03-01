function x = backsubst(A, b)
% BACKSUBST - Function that takes an upper triangular matrix A and a vector
%             of matching size b corresponding to the linear system Ax=b, 
%             and returns the solution for x.
%
%	MINIMAL WORKING EXAMPLE: Solve Ax = b for
%
%   A = [1 3 6; 0 5 2; 0 0 9];
%   b = [2.3; 1.4; 6.1];
%
%   x = backsubst(A,b); Solve Ax = b for x

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

%Sends error message if the given matrix is non-square
if size(A,1) ~= size(A,2)
   error('Matrix must be square'); 
end
%Sends error message if A and b are of incompatible sizes
if size(A,2) ~= size(b,1)
   error('Matrix/vector mismatch'); 
end

x = zeros(length(b),1); %Pre-define x as zero vector for speed
n = length(A); %Define number of rows in the square matrix A

%Starts at the last row and iteratively performs back substitution until it
%reaches the first row of the matrix. 
for i = n:-1:1
    x(i) = (b(i)-A(i,:)*x)/A(i,i);
end

%% Råkade göra kod för att transformera till upper triangular
% for i = 2:length(A)
% A(i:end,:) = A(i:end,:)-(A(i:end,i-1)/A(i-1,i-1))*A(i-1,:)
% end

