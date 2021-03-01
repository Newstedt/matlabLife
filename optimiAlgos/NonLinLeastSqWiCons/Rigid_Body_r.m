function [r,J,JJ]=Rigid_Body_r(x,p,w)
%RIGID_BODY_R - Function that setup residual function and its Jacobian for 
%               the photogrammetry optimization problem.
%
%   Minimum working example:
%
%   Used inside sqpsq:
%   [r,J,JJ]=Rigid_Body_r(x,p,w)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

%Perform cholesky factorization for handling weights
Rchol = chol(w);

%Unpack x to R
R={};
for i=0:2
    R{i+1} = [x((1:3)+i*9),x((4:6)+i*9),x((7:9)+i*9)]';
end

%Unpack x to d
d={};
for i=0:2
    d{i+1} = x((28:28+2)+i*3);
end

m = length(p)/3;

%Unpack x to q
k=0;
for j=1:m
       q(:,j) = [x(37+k*3);x(38+k*3);x(39+k*3)];
       k=k+1;
end

%Define residual function
for i=1:3
    r(1+3*m*(i-1):m*3+3*m*(i-1)) = reshape(R{i}*q + d{i} - p(:,(1:m)+m*(i-1)),[],1);
end

if nargout > 2 %We want the approximated Jacobian
    f = @(x) Rigid_Body_r(x,p);
    JJ = jacapprox(f,x);
    JJ = Rchol*JJ; %Take weights into consideration
end

%Compute the analytical Jacobian
J = zeros(length(r),length(x));
J1=[];
for i=1:m 
    J1 = [J1;blkdiag(q(:,i)',q(:,i)',q(:,i)')];
end
J1 = blkdiag(J1,J1,J1);
J2 = kron(ones(1,m)',eye(3));
J2 = blkdiag(J2,J2,J2);
J3 = [kron(eye(7)',R{1});kron(eye(7)',R{2});kron(eye(7)',R{3})];
J = [J1,J2,J3];
J = Rchol*J; %Take weights into consideration

r=r';
r = Rchol*r; %Take weights into consideration
