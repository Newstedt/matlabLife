np=2;
dim=3;
ns=150;
L=[1,0,0;1,0,0];
Ks=10;
Kd=0;
g=0;
Vini=[0, 0, 0; 0, 0, 0];
deltat=0.1;
m=1;
Xini=[0, 0, 0; 1.8, 0, 0];
V=Vini;

f=-Ks*(Xini-L);
x_old=Xini;
v_old=V-f*deltat/2;
x_new=x_old;
v_new=v_old;

triang=line([x_new(:,1);x_new(1,1);x_new(end,1)],[x_new(:,2);...
    x_new(1,2);x_new(end,2)],'marker','o','markersize',4);
set(triang,'erasemode','xor');

daspect([1,1,1]);
axis([-2,2,-2,2]);

hold on


for i=1:ns
    
    f=-Ks*(x_new-L);
    v_new=v_old+f*deltat/m;
    x_new=x_old+v_new*deltat;
    
    v_old=v_new;
    x_old=x_new;
    
    set(triang,'XData',[x_new(:,1);x_new(1,1);x_new(end,1)],...
        'YData',[x_new(:,2);x_new(1,2);x_new(end,2)]);
    pause(0.1);
end
