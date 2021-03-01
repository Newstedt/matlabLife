function [r,J,JJ] = circle_r_const(x,p)
%CIRCLE_R_CONST - Function that setup residual function for constrained 
%                 version of circle optimization problem.
%
%   Minimum working example:
%
%   Used inside sqpsq.m
%   [r,J] = circle_r_const(x,b)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

%Unpack x
p_i = x(4:end);

%construct residuals
p = reshape(p,[],1);
r = p_i - p;

if nargout>2 %We want the numerical jacobian
    f = @(x) circle_r_const(x,p);
    JJ = jacapprox(f,x);
end

if nargout>1 %We want the analytical jacobian
    J = [zeros(length(p),3), eye(length(p))];
end