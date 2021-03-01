function I = simpson(fun,M)
%MYQUAD approximates integral by using a composite quadrature rule
%
% I = MYQUAD(FUN,M) approximates the integral of scalar-valued
% function FUN from 0 to 1 by an appropriate quadrature rule that
% evaluates the function FUN in M points. FUN is a function handle
% and M is an odd positive integer. The function Y=FUN(X) should
% accept a vector argument X and return a vector result Y, the
% integrand evaluated at each element of X.

h=1/M;
x(1)=0;
I=0;

for i=1:M
    
    x(i+1)=x(i)+h;
    I=I+(h/6).*(fun(x(i))+4*fun((x(i)+x(i+1))/2)+fun(x(i+1)));

end


