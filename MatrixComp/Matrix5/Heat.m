function A=Heat(N,lambda)

% HEAT  Build matrix for the implicit Euler method for 2D heat equation
%
% INPUT:
%   N        the number of subintervals along each spatial dimension
%   lambda   the critical parameter linking k and h, k = lambda h^2
%
% OUTPUT:
%   A        a sparse matrix of dimension M^2, M = N-1 corresponding to 
%            Euler's method for the 2D heat equation
%
% MINIMAL WORKING EXAMPLE: Construct matrix of dimension 99*99
%
% A2=Heat(100,1);

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2017-10-10  Initial programming and testing
%   2018-10-22  Matrix B is now generated without wasting memory.

% Number of internal nodes per row
M=N-1;

% Sparse identity matrix of dimension M
I=speye(M,M);

% Construct the tridiagonal matrix B
%
%     B=lambda*toeplitz([2 -1 zeros(1,M-2)]); 
%
% as a sparse, rather than a dense matrix.

% Initalize arrrays to store COO representation of matrix B
i=zeros(3*M-2,1); j=zeros(3*M-2,1); val=zeros(3*M-2,1);

% Initialize counter to nonzero entries
k=1;

% Loop over rows
for t=1:M
    if t>1
        %Sub-diagonal entry
       i(k)=t; j(k)=t-1; val(k)=-1; k=k+1;
    end
    % Diagonal entry
    i(k)=t; j(k)=t; val(k)=2; k=k+1;
    if t<M
        %Super-diagonal entry
       i(k)=t; j(k)=t+1; val(k)=-1; k=k+1;
    end
end
% Scale the nonzero values.
val=val*lambda;

% Finally, generate the matrix B as a sparse matrix for i, j, val.
B=sparse(i, j, val);

% Construct the final matrix
A=kron(I,I)+kron(I,B)+kron(B,I);