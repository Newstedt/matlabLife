function [u] = g(x,y)

u=100*(0.25<x).*(x<0.75).*(0.25<y).*(y<0.75);
u=u+200*(0.1<x).*(x<0.20).*(0.45<y).*(y<0.55);
end
