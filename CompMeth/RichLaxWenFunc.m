function u = RichLaxWenFunc(N,M,fu,fx,x_start,x_end,lambda)

dx = (x_end-x_start)/N;
dt = lambda*dx;
x = dx/2:dx:x_end-dx/2;

u=zeros(N,M);

u(:,1)=fx(x);

plot(x,u(:,1));

for n = 1:M
    
        u_minus=2*1-u(1,n);
        u_plus=2*u(N,n)-u(N-1,n);
        u_new=[u_minus;u(:,n);u_plus];
    
        u_half=(1/2)*(u_new(1:N+1)+u_new([2:N+2]))-(1/2)*(dt/dx)*...
                    (fu(u_new([2:N+2]))-fu(u_new(1:N+1)));
        
        u(:,n+1)=u(:,n)-(dt/dx)*(fu(u_half([2:end]))-fu(u_half([1:end-1])));
        
        plot(x,u(:,n+1));
        axis([0 1 -1 3]);
        pause(.1);
end