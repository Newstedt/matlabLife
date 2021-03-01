function u = RichLaxWen(N,M,f,x_start,x_end,lambda,c)

dx = (x_end-x_start)/N;
dt = lambda*dx;
x = dx/2:dx:x_end-dx/2;

u=zeros(N,M);

u(:,1)=f(x);

plot(x,u(:,1));

for n = 1:M
        u_half(:,n)=(1/2)*(u(:,n)+u([2:N,1],n))-(1/2)*(dt/dx)*...
                    (c*u([2:N,1],n)-c*u(:,n));
        
        u(:,n+1)=u(:,n)-(dt/dx)*(c*u_half(:,n)-c*u_half([N,1:N-1],n));
        
        plot(x,u(:,n+1));
        axis([0 1 -1 1]);
        pause(.1);
end
