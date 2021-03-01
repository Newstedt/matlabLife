%Number of points, dimension, number of time steps and delta t
NP=2; dim=3; Nts=150; dt=0.1;

%Properties of string
NS=1; L=1; Ks=10; Kd=0;

%initial values
V0=[0, 0, 0; 0, 0, 0]; X0=[0, 0, 0; 1.8, 0, 0];
g=0; x_new=X0; v_new=V0; m=1;

spring_num=0;

%place the strings
for i=1:NS
    part_num=i;
    spring_num=spring_num+1;
    spring(spring_num).from=part_num;
    spring(spring_num).to=part_num+1;
    spring(spring_num).length=L;
    spring(spring_num).KS=Ks;
    spring(spring_num).KD=Kd;
end

f=zeros(NP,dim);
r=zeros(NP-1,dim);

for k=1:NP-1
      
        r(k,:)=x_new(k,:)-x_new(k+1,:);
        f(k,:)=-spring(k).KS*(norm(r(k,:))-L)*(r(k,:)/norm(r(k,:)));
        f(k+1,:)=-f(k,:);
    
end
% r=x_new(1,:)-x_new(2,:);
% f(1,:)=-spring(1).KS*(norm(r)-L)*(r/norm(r));
% f(2,:)=-f(1,:);

% f=-spring(1).KS*(norm(x_new)-L)*(x_new/norm(x_new));
    
x_old=x_new;
v_old=v_new-f*dt/2;

triang=line([x_new(:,1);x_new(1,1);x_new(end,1)],[x_new(:,2);...
    x_new(1,2);x_new(end,2)],'marker','o','markersize',4);
set(triang,'erasemode','xor');

daspect([1,1,1]);
axis([-0.1,1.9,-0.5,0.5]);
    
for i=1:Nts
    
%     f=-spring(1).KS*(norm(x_new)-L)*(x_new/norm(x_new));
    f=zeros(NP,dim);
    
    for k=1:NP-1
        
        r=zeros(NP-1,dim);
        r(k,:)=x_new(k,:)-x_new(k+1,:);
        f(k,:)=-spring(k).KS*(norm(r(k,:))-L)*(r(k,:)/norm(r(k,:)));
        f(k+1,:)=-f(k,:);
    
    end
        
    v_new=v_old+f*dt/m;
    x_new=x_old+v_new*dt;
    
    v_old=v_new;
    x_old=x_new;
    
    set(triang,'XData',[x_new(:,1);x_new(1,1);x_new(end,1)],...
        'YData',[x_new(:,2);x_new(1,2);x_new(end,2)]);
    pause(0.1);
end    
    
% for a=1:NP
%     
%     part_num=a;
%         X(part_num,1)=a*L;
%         X(part_num,2)=0;
%         X(part_num,3)=0;
%         V(part_num,1)=V0(1);
%         V(part_num,2)=0;
%         V(part_num,3)=0;
% end
