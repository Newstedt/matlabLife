function u = laxfriedpractice(g,f,x_0,x_end,M,N,lambda)

dx=(x_end-x_0)/M;
x=[x_0+dx/2:dx:x_end-dx/2];

u=g(x);
plot(x,u);

for n=1:N
    
    u = 1/2*(u([end,1:end-1])+u([2:end,1]))-1/2*lambda*(f(u([2:end,1]))-f(u([end,1:end-1])));
    
    plot(x,u);
    axis([0 1 -1 1]);
    pause(0.01);
    
end
