function alpha =linesearch(f,xk,pk,amin,F0,Fp0,c1)
%LINESEARCH Function that performs Armijo linesearch
%
%   Minimum working example:
%
%   Used inside gaussn.m
%   alpha =linesearch(f,xk,pk,amin,F0,Fp0,c1)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

j=0;
alpha = 2^-j;

% loop as long as alpha is larger than minimum allowed alpha
while(alpha > amin)
    % check convergence condition
    if (f(xk+alpha*pk) <= F0 + c1*alpha*Fp0)
        break
    end
    j=j+1; % step
    alpha = 2^-j; % bisect alpha for next iter
end