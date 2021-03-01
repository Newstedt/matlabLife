function u = laxfried(N,M,f,x_start,x_end,lambda,c)

dx = (x_end-x_start)/N;
dt = lambda*dx;
x = dx/2:dx:x_end-dx/2;

u=zeros(N, M);

u(:,1)=f(x);

plot(x,u(:,1));

for n = 1:M
        u(:,n+1)=(1/2)*(u([N,1:N-1],n)+u([2:N,1],n))-c.*(dt/(2*dx)).*(u([2:N,1],n)-u([N,1:N-1],n));
        
        plot(x,u(:,n+1));
        axis([0 1 -1 1]);
        pause(.1);
end
