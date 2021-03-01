%TEST_QR - Tests qr_fac.m
%
%   MINIMAL WORKING EXAMPLE: ">> test_qr" to see results in relevant tests.
%
% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

%Construct 4 different matrices of variable size
A1 = [-1 -1 1; 1 3 3; -1 -1 5; 1 3 7];
A2 = [1 1 1; 1 2 4; 1 3 9; 1 4 16];
A3 = rand(100,120);
A4 = rand(130,100);

%Call qr_fac for each of the matrices
[Q1,R1] = qr_fac(A1);
[Q2,R2] = qr_fac(A2);
[Q3,R3] = qr_fac(A3);
[Q4,R4] = qr_fac(A4);

%Check the relative error of the QR-factorisation for each case
rel_err(1) = norm(A1-Q1*R1)/norm(A1);
rel_err(2) = norm(A2-Q2*R2)/norm(A2);
rel_err(3) = norm(A3-Q3*R3)/norm(A3);
rel_err(4) = norm(A4-Q4*R4)/norm(A4);
rel_err = rel_err';

