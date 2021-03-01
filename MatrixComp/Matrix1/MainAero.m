% MAINAERO - Test script for 5 different scenarios of two air crafts leaving each known
%            location with a known speed.
%
%   MINIMAL WORKING EXAMPLE: 
%   ">> MainAero" to see results in relevant tests.

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

clear all; clc;
%% Define 5 different scenarios for p1, p2, v1, v2
p1_1 = [0.16; 0.02]; p2_1 = [0.35; 0];
v1_1 = [0.7; 0.84]; v2_1 = [-0.56; 0.79];

p1_2 = [0; 0]; p2_2 = [0.2; 0];
v1_2 = [1.6; 0.85]; v2_2 = [-1.78; 1.2];

p1_3 = [0; 0]; p2_3 = [1; 0];
v1_3 = [-1; 1]; v2_3 = [1; 1];

p1_4 = [0.01; 0]; p2_4 = [1.2; 0];
v1_4 = [1; 1]; v2_4 = [-1; 0.5];

p1_5 = [0; 0]; p2_5 = [0.35; 0];
v1_5 = [1; 1]; v2_5 = [-1.2; 2];

%% Call for collision function to find intersection point coordinates and
  %and at which time each air craft will be at the intersection point. 
  %This is done for each of the 5 scenarios.
[x1, y1, t1_1, t2_1] = collision(p1_1, v1_1, p2_1, v2_1);
[x2, y2, t1_2, t2_2] = collision(p1_2, v1_2, p2_2, v2_2);
[x3, y3, t1_3, t2_3] = collision(p1_3, v1_3, p2_3, v2_3);
[x4, y4, t1_4, t2_4] = collision(p1_4, v1_4, p2_4, v2_4);
[x5, y5, t1_5, t2_5] = collision(p1_5, v1_5, p2_5, v2_5);

%% Find longest time to intersetion point. This is for plotting purposes
tb_1 = max(t1_1,t2_1);
tb_2 = max(t1_2,t2_2);
tb_3 = max(t1_3,t2_3);
tb_4 = max(t1_4,t2_4);
tb_5 = max(t1_5,t2_5);

%% Cal trajCol.m for plotting the trajectories of the air crafts
trajCol(p1_1, v1_1, p2_1, v2_1, tb_1, 1);
trajCol(p1_2, v1_2, p2_2, v2_2, tb_2, 1);
trajCol(p1_3, v1_3, p2_3, v2_3, tb_3, 1);
trajCol(p1_4, v1_4, p2_4, v2_4, tb_4, 1);
trajCol(p1_5, v1_5, p2_5, v2_5, tb_5, 1);

%% Check which "cautious" level each scenario corresponds to
checkLevel(x1, y1, t1_1, t2_1);
checkLevel(x2, y2, t1_2, t2_2);
checkLevel(x3, y3, t1_3, t2_3);
checkLevel(x4, y4, t1_4, t2_4);
checkLevel(x5, y5, t1_5, t2_5);


