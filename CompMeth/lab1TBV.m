f=@(x) sin(2*pi*x);
u=laxfried(100,10,f,0,1,(49/100),2);

plot(u)

u_ans=RichLaxWen(100,10,f,0,1,(49/100),2);

figure
plot(u_ans)
