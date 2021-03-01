function [a] = newtonmet(r,p,x0,h,tol)

a=x0;
[f j]=soilf_J(a,r,p,h);

i=0;
    while norm(f)>tol && i<100
        i=i+1;
        s=-j\f';
        a=a+s.';
        [f j]=soilf_J(a,r,p,h);
    end
iter=i