function [Q,R] = qr_fac(A)
%QR_FAC - Takes matrix A and performs QR-factorisation on it. Returns Q & R
%
%   MINIMAL WORKING EXAMPLE: Perform QR-factorisation on a matrix A.
%   
%   A = [-1 -1 1; 1 3 3; -1 -1 5; 1 3 7]; %define A
%   [Q,R] = qr_fac(A); %perform QR-factorisation on A

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

m = size(A,1); %define the size of A
Q = eye(m); %pre-allocate Q
N = min(size(A)); %define how far the for-loop should loop

for i = 1:N
    
    %call gallery which returns components for householder matrix
    [v,beta] = gallery('house', A(i:end,i));
    %construct householder matrix from components received above
    H = eye(size(A(i:end,i:end),1))-beta*v*v';
    %overwrite A with householder matrix times A
    A(i:end,i:end) = H*A(i:end,i:end);
    
    %construct QA as an identity matrix of size m
    QA = eye(m);
    %set elements below current index element equal to householder equation
    QA(i:end,i:end) = H;
    %multiply Q with QA. This is to produce Q = Q_n*...*Q_2*Q_1*A.
    Q = Q*QA;
    
end
%set R equal to the remaining A
R = A;

