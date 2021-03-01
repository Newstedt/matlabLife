function [alpha] = linesearch(f,xk,pk,aMin,F0,Fp0,c1)

alpha = 1;

while aMin < alpha
    if feval(f,xk + alpha*pk) <= F0 + c1*alpha*Fp0
        break
    end
    alpha = alpha/2;
end