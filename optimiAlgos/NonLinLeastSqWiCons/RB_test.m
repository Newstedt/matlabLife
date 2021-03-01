%RB_TEST - Script that runs relevant tests for questions in
%             photogrammetry tasks in optimization assignment.
%
%   Minimum working example:
%
%   Used inside sqpsq:
%   >> RB_test

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.
clc
clear all
close all

%Set constants
maxIter=50; tol=1e-8; mu = 0.1; c1=0.1;aMin=1e-3;
stepType=0; nu0=0.1; epsC = 1e-8; epsR=1e-3;


%% Constraint problem
%Get points and information for MIT-building
for i=1:3
    [p,cTemp,cams]=mitpts(i);
    P(:,(1:7)+7*(i-1)) =p(:,1:7);
    c{i} = full(cTemp(1:21,1:21));
end

%Define weights
w1 = blkdiag(inv(c{1}),inv(c{2}),inv(c{3}));
w2 = eye(63);
w3 = blkdiag(eye(21),inv(c{2}),inv(c{3}));

%Define initial guess
R = reshape([eye(3);eye(3);eye(3)],[],1);
d =  zeros(9,1);
q = reshape(P(:,1:7),[],1);
x0 = [R;d;q];

%Define functions for residual/Jacobian- and constraint/Jacobian functions
fun = @(y) Rigid_Body_r(y,P,w3);
con = @(y) camG(y,P);

%Perform optimization using sqpsq
[x,n,code,l,X,a,C,L,nus,r,J,A]=sqpsq(fun,con,x0,epsR,epsC,maxIter,{},nu0,mu);

%Unpack x
R = reshape(x(1:27),3,[])';
d = x(28:36);
q = reshape(x(37:end),3,[]);

%Construct plots
for i=0:2
    figure
    hold all
    P_loc = R((1:3)+3*i,:)*q+d((1:3)+i*3);
    plot3(P_loc(1,:),P_loc(2,:),P_loc(3,:),'rx')
    plot3(P(1,(1:7)+i*7),P(2,(1:7)+i*7),P(3,(1:7)+i*7),'b*')
    title(sprintf("Camera: %d",i+1))
end