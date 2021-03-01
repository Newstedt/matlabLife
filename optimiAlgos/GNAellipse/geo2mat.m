function [U,b,gamma,Q]=geo2mat(c,ax,delta)
%GEO2MAT Convert ellipse equation parameters from geometrical to matrix form.
%
%[U,b,gamma,Q]=geo2mat(c,ax,delta)
%
%Input form:
%            (x-c(1))^2   (y-c(2))^2
%            ---------- + ---------- - 1 = 0
%              ax(1)^2      ax(2)^2
%with rotation delta (not shown).
%
%Output form:
%
% x'*U*x + b'*x - gamma = 0
%
%Q is the rotation matrix.

% v1.0  2000-05-13. Niclas Borlin, niclas@cs.umu.se.
% v1.1  2000-05-18. Added return of the rotation matrix.

Q=[cos(delta),-sin(delta);sin(delta),cos(delta)];
L=diag(1./ax.^2);

U=Q*L*Q';
b=-2*U*c;

gamma=1-c'*U*c;