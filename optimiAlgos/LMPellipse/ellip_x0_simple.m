function [c,ax,delta,theta] = ellip_x0_simple(x,y)

c(1,1) =   mean(x);
c(2,1) =   mean(y);
ax = [abs((max(x)-min(x)))/2;abs((max(y)-min(y)))/2];
theta = linspace(0,2*pi,length(x));
delta = 0;