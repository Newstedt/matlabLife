function [I] = simpson2(f,h,xspan)

I=0;
n=(xspan(end)+xspan(1))/h;

for i=1:3
    
    I=I+(h/3).*(f(2*i-1)+4.*f(2*i)+f(2*i+1));
  
end

