function[] = checkLevel(x, y, t1, t2)
%CHECKLEVEL - Function that checks alert Levels 1 - 5 for two air plane 
%             trajectories, and gives corresponding error messages.
%
%   Input parameters:
%   x: x-coordinate for intersection point of two air craft trajectories
%   y: y-coordinate for intersection point of two air craft trajectories
%   t1: Time at which the fastest air craft reaches the intersection point
%   t2: Time at which the slowest air craft reaches the intersection point
%
%   MINIMAL WORKING EXAMPLE: Check levels for two airplanes with
%   intersection point at x = 2.3 and y = 4.8, and times of intersection
%   t1 = 8.92 and t2 = 38.2.
%
%   >> checkLevel(x, y, t1, t2)
%   >> "Gives level of alert"

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-13: Initial version 
% 2018-09-27: Edited comments
%
% Function code starts here...

fprintf('Check result: \n');
    
dt = abs(t1-t2); %Compute dt = time difference between t1 and t2 
dtM = floor(dt*60); %Convert time difference between t1 and t2 to minutes
dtS = (dt*60-dtM); %Seconds after last minute in time difference
ta = min(t1,t2);
taM = floor(ta*60);
taS = (ta*60-taM)*60;

%Level 1 - start by checking if the paths are even crossing each other in 
%the future. If not, all other levels are irelevant.
if isnan(x) || isnan(y)
    fprintf('NO CONFLICT POSSIBLE \n \n');

%Check condition for each level

%Level 5
elseif ta <= 3*(1/60) && dt <= 3*(1/60)
    
    fprintf('CONFLICT! PREDICTED COLISION IN %gM%02gS (%gM%02gS) \n \n',...
        taM, taS, dtM, dtS);
    
%Level 4
elseif ta <= 1/6 && dt <= 3*(1/60)
    
    fprintf('CONFLICT IN %gM%02gS (%gM%02gS) \n \n', taM-3, taS,...
        dtM, dtS);
    
%Level 3
elseif ta <= 1/6 && dt > 3*(1/60)

    fprintf('NEAR CONFLICT IN %gM%02gS (%gM%02gS) \n \n', taM, taS,...
            dtM, dtS);
    
%Level 2
elseif ta > 1/6 || dt > 1/6
    fprintf('NO CONFLICT \n \n');    
end
end