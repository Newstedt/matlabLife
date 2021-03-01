function u = laxfriedVectEdit(N,M,fu,fx,x_start,x_end,lambda,base,increase)

dx = (x_end-x_start)/N;
dt = lambda*dx;
x = dx/2:dx:x_end-dx/2;

u=[fx{1}(x); fx{2}(x)];

plot(x,u(1,:));

h=base+increase*(x<=120);

for n = 1:M
    
    u1=(1/2).*([u(1,2); u(2,2)]+[u(1,1); u(2,1)])...
        -(dt/(2*dx))*...
        ([fu{1}((u(1,2)),(u(2,2)));...
          fu{2}((u(1,2)),(u(2,2)))]...
        -[fu{1}((-u(1,1)),(-u(2,1)));
          fu{2}((-u(1,1)),(-u(2,1)))]);
    
    uN=(1/2).*([u(1,N); u(2,N)]+[u(1,N-1); u(2,N-1)])...
        -(dt/(2*dx))*...
        ([fu{1}((-u(1,N)),(-u(2,N)));...
          fu{2}((-u(1,N)),(-u(2,N)))]...
        -[fu{1}((u(1,N-1)),(u(2,N-1)));
          fu{2}((u(1,N-1)),(u(2,N-1)))]);
    
        uMid=(1/2)*([u(1,[1:N-2]); u(2,[1:N-2])]+[u(1,[3:N]); u(2,[3:N])])...
           -(dt/(2*dx)).*...
           ([fu{1}(([u(1,[3:N])]),([u(2,[3:N])]));...
             fu{2}(([u(1,[3:N])]),([u(2,[3:N])]))]...
           -[fu{1}(([u(1,[1:N-2])]),([u(2,[1:N-2])]));
             fu{2}(([u(1,[1:N-2])]),([u(2,[1:N-2])]))]);
         
         u=[u1,uMid,uN];
         
%          h=base+increase*(x<=120);

         position(n)=dx*find(u(1,:)>base+increase/4,1,'last');
         
         
        plot(x,u(1,:));
        axis([0 200 -1 50]);
        pause(.0001);
end

v=(position(202)-position(2))/(200*dt)

