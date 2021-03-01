%TEST Script which runs relevant tests for Assignment 2 in the course  
%     "Engineering Optimization" by calling test_optimization.
%
%   Minimum working example:
%   
%   >> Test

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

clc
clear all
close all
%% ELLIPSE PROBLEM Gauss, Naive x0, [1,1], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,1],'Gauss','Naive');
[r,J] = Residual_Function(x,p);
sigma0 = sqrt(r'*r/(length(p)))
c_xx = sigma0^2*inv(J'*J)
%% ELLIPSE PROBLEM Gauss, Naive x0, [1,1], failure
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,2],'Gauss','Naive');

%%
%% ELLIPSE PROBLEM Gauss, Ellip x0, [1,1], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,1],'Gauss','ellip_x0');

%% ELLIPSE PROBLEM Gauss, Ellip x0, [1,2], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,2],'Gauss','ellip_x0');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELLIPSE PROBLEM LM, Naive x0, [1,1], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([2,2],'LM','Naive');

%% ELLIPSE PROBLEM LM, Naive x0, [1,1], failure
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,2],'LM','Naive');

%% ELLIPSE PROBLEM LM, Naive x0, [1,1], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([2,2],'LM','ellip_x0');

%% ELLIPSE PROBLEM LM, Naive x0, [1,1], completed
[x,code,n,X,dk,rhos,stepType,p] = test_optimization([1,2],'LM','ellip_x0');
