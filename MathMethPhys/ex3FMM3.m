%Number of points, dimension, number of time steps and delta t
NP=32; dim=3; Nts=1500; dt=0.01; Num_col=8; Num_row=4;

%Properties of string
NS=1; L=1; Ks=10; Kd=0;

%initial values
V0=[0, -5, 0; 0, 5, 0]; X0=[0, 0, 0; 1.8, 0, 0];
g=[0, 0, 0; 0, 0, 0]; x_new=X0; v_new=V0; mval=[1 1]; m=diag(mval);

spring_num=0;

for a=1:NP
    
    part_num=a;
        X(part_num,1)=a*L;
        X(part_num,2)=a-1;
        X(part_num,3)=a;
        X(part_num,4)=0;
        V(part_num,1)=V0(1);
        V(part_num,2)=0;
        V(part_num,3)=0;
end

%place the strings
    spring_num=spring_num+1;
    for j=1:Num_col-1
        part_num=i;
        spring(j).from=part_num;
        spring(j).to=part_num+1;
        spring(j).length=L;
        spring(j).KS=Ks;
        spring(j).KD=Kd;
        
end

f=zeros(NP,dim);
r=zeros(NP-1,dim);

for k=1:NP-1
      
        r(k,:)=x_new(k,:)-x_new(k+1,:);
        rdot=v_new(k,:)-v_new(k+1,:);
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
axis([-0.5,2.9,-1.5,1.5]);
    
for i=1:Nts
    
%     f=-spring(1).KS*(norm(x_new)-L)*(x_new/norm(x_new));
    f=zeros(NP,dim);
    
    for k=1:NP-1
        
        r=zeros(NP-1,dim);
        rdot(k,:)=v_new(k,:)-v_new(k+1,:);
        r(k,:)=x_new(k,:)-x_new(k+1,:);
        f(k,:)=-(spring(k).KS*(norm(r(k,:))-L)+...
            spring(k).KD*((sum(conj(rdot(k,:)).*r(k,:)))/norm(r(k,:))))...
            *(r(k,:)/norm(r(k,:)));
        f(k+1,:)=-f(k,:);
        fjf(k)=norm(r(k,:))-L;
        amp(i)=fjf(k);
    end
      
    v_new=v_old+m\f*dt;
    x_new=x_old+v_new*dt;
    v_k=(v_old+v_new)/2; %mittenvärde
    
    v_old=v_new;
    x_old=x_new;
    
 Ek=sum(sum((1/2)*(m*(v_k).^2)));
 Ep=sum((m*g)'*x_new(:,2));
 Efj=(1/2)*spring.KS*fjf.^2;
 Etot(i)=Ek+Ep+Efj;
 
 dist_cent=x_new-0.9;
 ang_moment(i)=sum(m*(dist_cent(:,1).*v_new(:,2)-dist_cent(:,2).*v_new(:,1)));
 
 W(i)=sum(sqrt((v_new(:,1)).^2+(v_new(:,2)).^2));
 spring_length(i)=norm(r);
 
    set(triang,'XData',[x_new(:,1);x_new(1,1);x_new(end,1)],...
        'YData',[x_new(:,2);x_new(1,2);x_new(end,2)]);
%      pause(0.01);
end    

diff=max(Etot)/min(Etot);
figure
plot(Etot)
figure
plot(amp)

j=0;

if Kd>0
    
for i=1:Nts-1
    
    amp2(i)=abs(amp(i));
    
    if(amp(i)<0 && amp(i+1)>0 || amp(i)>0 && amp(i+1)<0)
       
        j=j+1;
        maxamp(j)=max(amp2);
        amp2=0;
        
        if maxamp(j)<0.1*maxamp(1)
            percent_limit=i
            break
        end
    end
end
plot(maxamp);
end

figure
ang_moment(end);
plot(ang_moment)

figure
plot(W)

figure
plot(spring_length)
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