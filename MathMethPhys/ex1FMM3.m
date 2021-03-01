np=3;
dim=3;
ns=1500;
L=1;
Vini=0.01;
deltat=0.1;
m=1;
X=[-L/2, 0, 0 ; L/2, 0, 0 ; 0, L/2, 0];
V=[Vini, 0, 0 ; 0, -Vini, 0 ; -Vini/2, -Vini/2, 0];

f=0.5*(0.5-rand(np,dim))-0.1*V;
x_old=X;
v_old=V-f*deltat/2;
x_new=x_old;
v_new=v_old;

triang=line([x_new(:,1);x_new(1,1);x_new(end,1)],[x_new(:,2);...
    x_new(1,2);x_new(end,2)],'marker','o','markersize',4);
set(triang,'erasemode','xor');

daspect([1,1,1]);
axis([-5,5,-5,5]);

hold on


for i=1:ns
    
    f=0.5*(0.5-rand(np,dim))-0.1*v_new;
    v_new=v_old+f*deltat/m;
    x_new=x_old+v_new*deltat;
    
    v_old=v_new;
    x_old=x_new;
    
    set(triang,'XData',[x_new(:,1);x_new(1,1);x_new(end,1)],...
        'YData',[x_new(:,2);x_new(1,2);x_new(end,2)]);
    pause(0.01);
end


