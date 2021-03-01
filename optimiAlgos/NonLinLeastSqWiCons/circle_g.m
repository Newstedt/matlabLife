function pts = circle_g(c,r,theta)
%CIRCLE_G - Function that setup points from model function for circle 
%           optimization problem.
%
%   Minimum working example:
%
%   Used inside circle_r:
%   pts = circle_g(c,r,th)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

% Adjust theta and get length
theta=theta(:)';
n = length(theta);

% Construct points from model function
pts=zeros(2,n);
for i=1:n
    pts(:,i)=c+r*[cos(theta(i));sin(theta(i))];
end