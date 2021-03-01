function [c,A,AA] = circle_gx_const(x,p)
%CIRCLE_GX_CONST - Function that setup constraint function  and its 
%                  Jacobian for constrained version of the circle 
%                  optimization problem.
%
%   Minimum working example:
%
%   Used inside sqpsq.m:
%   [c,A,AA] = circle_gx_const(x,p)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

%Unpack x
c = x(1:2);
r = x(3);
p_i = reshape(x(4:end),[],length(x(4:end))/2);

%Define constraint function
for i = 1:length(p_i)
    gx(i) = (p_i(:,i) - c)'*(p_i(:,i) - c) - r^2;
end
    
if nargout>2 %We want the numerical jacobian
    f = @(x) circle_gx_const(x,p);
    AA = jacapprox(f,x);
end

A = zeros(length(p),2*length(p)+3);

if nargout>1 %We want the analytical jacobian
    for i = 1:length(p)
        A(i,:) = [-2*(p_i(:,i)-c)', -2*r, zeros(1,2*length(p_i))];
        A(i,(1:2)+2*(i-1)+3) = 2*(p_i(:,i)-c)';
    end
end

c = gx'; %Flip gx and set to c for return