function u = laxfriedpractice2(g,f,x_0,x_end,M,N,lambda,base,increase)

dx=(x_end-x_0)/M;
x=[x_0+dx/2:dx:x_end-dx/2];

u=g(x);

plot(x,u(1,:));

for n=1:N
    
    u_first = 1/2*(u(:,2)+u(:,1))-1/2*lambda*(f(u(:,2))-f(-u(:,1)));
    u_last = 1/2*(u(:,end)+u(:,end-1))-1/2*lambda*(f(-u(:,end))-f(u(:,end-1)));
    
    u=[u_first,u,u_last];
    
    u = 1/2*(u(:,1:end-2)+u(:,3:end))-1/2*lambda*(f(u(:,3:end))-f(u(:,1:end-2)));
    
    position(n) = dx*find(u(1,:)>base+increase/4,1,'last');
    
    plot(x,u(1,:));
    axis([0 200 0 40]);
    pause(0.0001);
    
end

velocity=(position(202)-position(2))/(200*lambda*dx)
