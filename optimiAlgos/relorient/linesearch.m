function alpha =linesearch(f,xk,pk,amin,F0,Fp0,c1)
j=0;
alpha = 2^-j;

while(alpha > amin)
    if (f(xk+alpha*pk) <= F0 + c1*alpha*Fp0)
        break
    end
    j=j+1;
    alpha = 2^-j;
end