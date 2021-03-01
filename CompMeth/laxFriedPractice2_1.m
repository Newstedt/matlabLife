function u = laxFriedPractice2_1(M,N,x0,xend,f,g,lambda)

dx=(xend-x0)/M;
x=[x0+dx/2:dx:xend-dx/2];

u=g(x);

plot(x,u(1,:));

for n=1:N
    
    u_first=1/2*(u(:,2)+u(:,1))-1/2*lambda*(f(u(:,2))-f(-u(:,1)));
    u_last=1/2*(u(:,end)+u(:,end-1))-1/2*lambda*(f(-u(:,end))-f(u(:,end-1)));
    
    u=[u_first,u,u_last];
    
    u=1/2*(u(:,1:end-2)+u(:,3:end))-1/2*lambda*(f(u(:,3:end))-f(u(:,1:end-2)));
    
end
    
    
    
