g=9.82;
base=15;
increase=15;
hq0={@(x) (base+increase).*((x<=120)) + base.*((120<x)); @(x) 0*x};

fhq= {@(h,q) q; @(h,q) (q.^2)./h+g*(h.^2)/2};

u = LaxFriedVectEdit(100,1000,fhq,hq0,0,200,0.05,base,increase);

% x=linspace(0,200);
% 
% plot(x,hq0{1}(x));

