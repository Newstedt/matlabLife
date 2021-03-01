function [r,J] = Residual_Function(x,p)
%Residual_Function - Function that setup residual function for optimization
%                    problems
%
%   Minimum working example:
%
%   Used inside gaussn.m
%   [r,J] = Residual_Function(x,p)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

h = Model_Function(x); % get model function
r = h-p; % compute residuals
r = reshape(r,[],1); % construct r as column vector
J = zeros(2*length(x(6:end)),length(x(6:end))+ 5); % pre-set Jacobi matrix

% compute Jacobi matrix
for i = 1:length(x(6:end))
    J(2*i-1:2*i,:) = ...
    [1 0 cos(x(5))*cos(x(5+i)) -sin(x(5))*sin(x(5+i)) ...
    -x(3)*sin(x(5))*cos(x(5+i))-x(4)*cos(x(5))*sin(x(5+i))...
    zeros(1,length(x(6:end)));...
    0 1 sin(x(5))*cos(x(5+i)) cos(x(5))*sin(x(5+i)) ...
    x(3)*cos(x(5))*cos(x(5+i))-x(4)*sin(x(5))*sin(x(5+i))...
    zeros(1,length(x(6:end)))];...
    J(2*i-1:2*i,i+5) = [-x(3)*cos(x(5))*sin(x(5+i))...
    -x(4)*sin(x(5))*cos(x(5+i));...
    -x(3)*sin(x(5))*sin(x(5+i))+x(4)*cos(x(5))*cos(x(5+i))];
end
end