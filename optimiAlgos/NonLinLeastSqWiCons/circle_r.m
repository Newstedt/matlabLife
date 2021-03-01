function [r,J,JJ] = circle_r(x,b)
%CIRCLE_R - Function that setup residual function for circle optimization
%           problem
%
%   Minimum working example:
%
%   Used inside gaussn.m
%   [r,J] = circle_r(x,b)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

%Unpack x
c = x(1:2);
r = x(3);
th = x(4:end);

%verify sizes
if length(th) ~= size(b,2), error('Wrong size'); end

%Let circle_g compute the points. Unroll difference to column vector.
r = reshape(circle_g(c,r,th)-b,[],1);

if nargout>2 %We want the numerical jacobian
    f = @(x) circle_r(x,b);
    JJ = jacapprox(f,x);
end

if nargout>1 %We want the analytical jacobian
    for i = 1:length(x(4:end))
        J(2*i-1:2*i,:) = [1 0 cos(x(3+i)) zeros(1,length(x(4:end)));
                          0 1 sin(x(3+i)) zeros(1,length(x(4:end)))];
        
        J(2*i-1:2*i,i+3) = [-x(3)*sin(x(3+i)); x(3)*cos(x(3+i))];
    end
end