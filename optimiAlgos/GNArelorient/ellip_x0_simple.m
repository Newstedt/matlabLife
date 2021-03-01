function [c,ax,delta,theta] = ellip_x0_simple(x,y)
%ELLIP_X0_SIMPLE - Function that computes a naive initial guess for the 
%                  ellipse problem.
%
%   Minimum working example:
%
%   Used inside test_optimization.m
%   [c,ax,delta,theta] = ellip_x0_simple(x,y)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

c(1,1) = mean(x);
c(2,1) = mean(y);
ax = [abs((max(x)-min(x)))/2;abs((max(y)-min(y)))/2];
theta = linspace(-pi,pi,length(x));
delta = 0;