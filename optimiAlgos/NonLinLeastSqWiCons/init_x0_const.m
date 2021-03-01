function x0 = init_x0_const(p)
%INIT_X0_CONST - Function that setup an initial guess for the constrained
%                circle optimization problem.
%
%   Minimum working example:
%
%   Call with problem specific 2D-points vector p:
%   x0 = init_x0_const(p)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

cx = mean(p(1,:)); %initial value for cx
cy = mean(p(2,:)); %initial value for cx

%initial value for r
r = mean([(max(p(1,:))-min(p(1,:))), (max(p(2,:))-min(p(2,:)))])/2;

%initial value for p_i
p_i = reshape(p,[],1);

%gather in one x0-vector for returning
x0 = [cx,cy,r,p_i']';



