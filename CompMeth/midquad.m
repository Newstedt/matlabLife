function I = midquad(fun,M)
%MIDQUAD approximates integral by using midpoint rule
%
% I = MIDQUAD(FUN,M) approximates the integral of the
% scalar-valued function FUN from 0 to 1 by using the composite
% midpoint rule on M subintervals. FUN is a function handle and
% M is a positive integer. The function Y=FUN(X) should accept
% a vector argument X and return a vector result Y, the
% integrand evaluated at each element of X. 
h=1/M;
t=0;
I=0;

for i = 1:M
    
    I = I + h.*fun(t + h/2);
    
    t=t+h;

end