function [x, y, t1, t2] = collision(p1, v1, p2, v2)
%COLLISION - Function that takes the initial coordinates and velocities of 
%            two air crafts and returns the intersection point of their 
%            trajectories as well as the time at which each air craft 
%            reaches this location.
%
%   Input parameters:
%   p1: x- & y-coordinates of air craft 1's initial position
%   v1: directional velocity of air craft 1
%   p2: x- & y-coordinates of air craft 2's initial position
%   v2: directional velocity of air craft 2
%
%   Output parameters:
%   x: x-coordinate of intersection point
%   y: y-coordinate of intersection point
%   t1: Time at which air craft 1 reaches the intersection point
%   t2: Time at which air craft 2 reaches the intersection point
%
%   MINIMAL WORKING EXAMPLE: Find intersection point coordinates and times
%   of intersection for two air crafts with initial coordinates:
%   p1 = [0.16; 0.02]; p2 = [0.35; 0];
%   and initial velocities: 
%   v1 = [0.7; 0.84]; v2 = [-0.56; 0.79];
%
%   [x, y, t1, t2] = collision(p1, v1, p2, v2)

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

%Construct matrix A with v1 and -v2 as column 1 and 2
A = [v1, -v2]; 
%Construct vector b from subtracting position coordinates p1 from position
%coordinates p2
b = p2-p1;
%Solve the system At=b using backslash operator
t = A\b;
%Give value to t1 and t2 so that the function can return them
t1 = t(1); t2 = t(2);

x = p1(1) + t(1)*v1(1); %Calculate x (for return)
y = p1(2) + t(1)*v1(2); %Calculate y (for return)

%Check if the air craft trajectories will cross in the future.
if t1 < 0 || t2 < 0
   x = NaN; y = NaN; %If not, return NaN for their coordinates 
end

end

