fx = @(x) 0*(x>0.5) + 1*(x<=0.5);

fu = @(u) 2*u;

u=laxfriedFunc(100,50,fu,fx,0,1,(49/100));

u_ans=RichLaxWenFunc(100,50,fu,fx,0,1,(49/100));
