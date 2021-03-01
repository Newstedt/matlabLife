function u = richlaxwenPracticeGhost2(M,N,x0,xend,f,g,lambda)

dx=(xend-x0)/M;
x=x0+dx/2:dx:xend-dx/2;

u=g(x);

plot(x,u);

for n=1:N
    
    u0 = 2-u(1);
    uMplus=2*u(end)-u(end-1);
    u=[u0,u,uMplus];
    
    u_half = (1/2)*(u(1:end-1)+u(2:end))-(1/2)*lambda*(f(u(2:end))-f(u(1:end-1)));
    
    u = u(2:end-1) - lambda*(f(u_half(2:end))-f(u_half(1:end-1)));
    
    plot(x,u);
    axis([0 1 -1 1]);
    pause(0.1);
    
end
    
    
    