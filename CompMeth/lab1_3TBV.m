fx= @(x) 1+0.7*sin(10*x);
fu= @(u) u.^2/2;

u_ans=RichLaxWenFunc(100,50,fu,fx,0,1,(49/100));
