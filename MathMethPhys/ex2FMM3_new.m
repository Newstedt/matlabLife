clear all
clc
%Number of points, dimension, number of time steps and delta t
NP=2; dim=3; Nts=700; dt=0.009;

%Properties of string
NS=1; L=1; Ks=10; Kd=0;

%initial values
V0=[0, 5, 0; 0, -5, 0]; X0=[0, 0, 0; 1.8, 0, 0];
g=[0, 0, 0; 0, 0, 0]; x_new=X0; v_new=V0; mval=[1 1]; m=diag(mval);
cent=(X0(2,:)-X0(1,:))/2; p=0;

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
        rdot=v_new(k,:)-v_new(k+1,:);
        f(k,:)=-(spring(k).KS*(norm(r(k,:))-L)+...
            spring(k).KD*((sum(conj(rdot(k,:)).*r(k,:)))/norm(r(k,:))))...
            *(r(k,:)/norm(r(k,:)));
        f(k+1,:)=-f(k,:);
    
end
    
x_old=x_new;
v_old=v_new-f*dt/2;

triang=line([x_new(:,1);x_new(1,1);x_new(end,1)],[x_new(:,2);...
    x_new(1,2);x_new(end,2)],'marker','o','markersize',4);
set(triang,'erasemode','xor');

daspect([1,1,1]);
axis([-0.95,2.9,-1.5,2]);
%------ make video file ------
% writerObj = VideoWriter('rotating_spring.avi');
% writerObj.FrameRate = 60;
% open(writerObj);
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
%---- continue after loop -----

% W=zeros(NP,dim); 

for i=1:Nts
    
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
    
 Ek(i)=sum(sum((1/2)*(m*(v_k).^2)));
 Ep(i)=sum((m*g)'*x_new(:,2));
 Efj(i)=(1/2)*spring.KS*fjf.^2;
 Etot(i)=Ek(i)+Ep(i)+Efj(i);
 
 dist_cent=0.9;
 ang_moment(i) = sum(m*(x_new(:,2).*v_new(:,1)-(x_new(:,1)-dist_cent).*v_new(:,2)));
 W(i)=ang_moment(i)/(fjf^2*m(1));
 
 spring_length(i)=norm(r);

    set(triang,'XData',[x_new(:,1);x_new(1,1);x_new(end,1)],...
        'YData',[x_new(:,2);x_new(1,2);x_new(end,2)]);
%      pause(0.01);

%---- part of movie making ----    
%     frame = getframe;
%    writeVideo(writerObj,frame);
%------------------------------
end 
%---- end movie making ------
% close(writerObj);
%----------------------------

% timesteps_per_lap=2*(sum(timesteps_half_lap)/p);
% time_per_lap=dt*timesteps_per_lap

diff=max(Etot)/min(Etot)
figure
plot(Etot)
hold on
plot(Ek)
hold on
plot(Ep)
hold on
plot(Efj)
legend('Etot','Ek','Ep','Efj')
% figure
% plot(amp)


if Kd>0

    for i=1:Nts-1;
    
        if amp(Nts-i) > 0.08;
            k=i;
            time2 = Nts * dt - k*dt;
            break
        end
    end
    time2
end

equi=sum(spring_length)/Nts

figure
ang_moment(end);
plot(ang_moment)
legend('angular momentum')

figure
plot(W)
hold on

plot(spring_length)

legend('angular freq','spring length')
