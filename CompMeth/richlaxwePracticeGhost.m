function u = richlaxwePracticeGhost(g,f,x_0,x_end,M,N,lambda)

dx=(x_end-x_0)/M;
x=[x_0+dx/2:dx:x_end-dx/2];

u=g(x);

plot(x,u);

for n = 1:N
    
    u_minus = 2-u(1);
    u_plus = 2*u(end)-u(end-1);
    
    u = [u_minus,u,u_plus];
    
    u_half = 1/2*(u([1:end-1])+u([2:end]))-1/2*lambda*(f(u([2:end]))-f(u([1:end-1])));
    
    u = u([2:end-1]) - lambda*(f(u_half(2:end))-f(u_half(1:end-1)));
    
    plot(x,u);
    axis([0 1 -1 2]);
    pause(0.01);
    
end