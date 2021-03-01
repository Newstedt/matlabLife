
function[] = trajCol(p1, v1, p2, v2, tb, print)
% TRAJCOL - Function that takes the initial positions and velocities of two
%           air crafts and plots their trajectories until they have both 
%           reached theirintersection point (if there is one, else it plots
%           one hour of flying).
% 
%   Iput parameter:
%   p1: x- & y-coordinates of air craft one's initial position
%   v1: directional velocity of air craft one
%   p2: x- & y-coordinates of air craft two's initial position
%   v2: directional velocity of air craft two
%   tb: time at which the slower air craft reaches the intersection point
%   print: set = 1 if you want to plot set = 0 if you don't want to plot
%
%   MINIMAL WORKING EXAMPLE: Plot trajectories of two air crafts with
%   initial positions:
%   p1 = [0.16; 0.02]; p2 = [0.35; 0];
%   initial velocities: 
%   v1 = [0.7; 0.84]; v2 = [-0.56; 0.79];
%   intersection time for second air craft:
%   tb = 17.4;
%   and print = 1.
%
%   >> trajCol(p1, v1, p2, v2, tb, print)
%   >> *plots trajectories*

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

t = 0:1/6000:tb; %Set time to slower plane reaching intersection point
tbM = floor(tb*60); %Convert tb to integer minutes
tbS = (tb*60-tbM)*60; %Convert tb to integer seconds over last minute

%In case of no intersection point in positive time
if tb < 0
    t = 0:1/6000:1;
end

%Compute trajectory for air craft one
x1 = p1(1) + v1(1)*t;
x1 = x1';
y1 = p1(2) + v1(2)*t;
y1 = y1';
pos1 = [x1 y1];

%Compute trajectory for air craft two
x2 = p2(1) + v2(1)*t;
x2 = x2';
y2 = p2(2) + v2(2)*t;
y2 = y2';
pos2 = [x2 y2];

%Check if the air crafts has intersection time in positive time and if user
%want to plot or not.
if tb > 0 && print == 1
    %Plot the trajectories
    figure
    hold on
    plot(pos1(:,1), pos1(:,2),'-.','LineWidth',2);
    plot(pos2(:,1), pos2(:,2),'-.','LineWidth',2);
    title(['Air crafts trajectories after ',num2str(tbM),'M and ',num2str(tbS),'S']);

%Check if the air crafts has intersection time in negative time and if user
%want to plot or not.
elseif tb < 0 && print == 1
    %Plot the trajectories
    figure
    hold on
    plot(pos1(:,1), pos1(:,2),'-.','LineWidth',2);
    plot(pos2(:,1), pos2(:,2),'-.','LineWidth',2);
    title('Air crafts trajectories that never cross each other');    
end

end