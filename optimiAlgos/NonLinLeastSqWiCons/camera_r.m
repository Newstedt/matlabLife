function [r,J,JJ] = camera_r(x,p)
%CAMERA_R - Function that setup residual function for optimization
%           problems
%
%   Minimum working example:
%
%   Used inside gaussn.m
%   [r,J] = Residual_Function(x,p)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

%Shortcut for selftest
if ischar(x), selftest, return; end

%Unpack x
R = x(1:27);
d = x(28:36);
q = reshape(x(37:end),3,[]);
r = zeros(3*length(q),1);

for i = 1:3
    R_i = reshape(R((i-1)*9+1:9*i),[],3)';
    U = R_i*q + d(1+(i-1)*3:3*i).*ones(3,length(q));
    qTil = p(:,(i-1)*length(p)/3+1:i*length(p)/3);
    r((i-1)*length(p)+1:i*length(p)) = reshape(U - qTil,[],1);
end

r = reshape(r,[],1);

if nargout>2 %We want the numerical jacobian
    f = @(x) camera_r(x,p);
    JJ = jacapprox(f,x);
end

if nargout>1 %We want the analytical jacobian
   J = 0;
end