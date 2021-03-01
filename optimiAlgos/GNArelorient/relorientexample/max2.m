function [m,r,c]=max2(A)
%MAX2 Returns maximum element of a matrix
%
%[m,r,c]=max2(A)
%m is the maximum value
%r is the row of the maximum value
%c is the column of the maximum value

% v1.0  94-05-07. Niclas Borlin, niclas@cs.umu.se.
% v1.1  94-12-12. Modified to also return position of maximum value.

[m,r]=max(A,[],1);
[m,c]=max(m);
r=r(c);
