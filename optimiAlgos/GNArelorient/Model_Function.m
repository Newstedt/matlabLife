function h = Model_Function(x)
%Model_Function - Function that setup model function for optimization
%
%   Minimum working example:
%
%   Used inside Residual_Function.m
%   h = Model_Function(x)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

% define rotation matrix 
Q = [cos(x(5)) -sin(x(5)); sin(x(5)) cos(x(5))];

% setup h for all points
for i = 1:length(x(6:end))
    h(i,:) = x(1:2) + Q*[x(3) 0; 0 x(4)]*[cos(x(5+i));sin(x(5+i))];
end
h = h';
end