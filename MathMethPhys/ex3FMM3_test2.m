clc
clear all
%Number of points, dimension, number of time steps and delta t
NP=32; dim=3; Nts=1000; dt=0.002; Num_col=8; Num_row=4; r=zeros(1,dim);
rdot=zeros(1,dim); f=zeros(NP,dim); gravity=1; g=zeros(NP,dim); 
g(:,2)=-gravity; 
% n=1;

%Properties of string
NS=96; L=1; L2=sqrt(2); Ks=100; Kd=0;

%initial values
mval=ones(NP,1); m=diag(mval); v0=[5,0,0];

%initial values balls
NumBalls=16;
X_ball=0;
Y_ball=0;
R_ball=L;

%position-vector first ball
Pos_ball=zeros(NumBalls,dim);
Pos_ball(1,:)=[X_ball, Y_ball, 0];

%placing balls
% +0.1*R_ball*randn
% +[R_ball+0.1*R_ball*randn,0,0]
for i=2:NumBalls
    Pos_ball(i,:)=Pos_ball(i-1,:)+[R_ball+0.1*R_ball*randn,0,0];
    BOLL_a=rectangle('Position',[Pos_ball(i,1)-R_ball,Pos_ball(i,2)-R_ball,2*R_ball,2*R_ball],...
    'Curvature',[1,1],'EdgeColor','r');
end

spring_num=0;
k=0;

for i=1:Num_row
    for j=1:Num_col
        k=k+1;
        x(k,:)=[j,i,0];
        v(k,:)=[5,0,0];
    end
end


x_new=x;
v_new=v;
x_old=x_new;

part_num=0;

%place the strings
for j=1:Num_row
    part_num=part_num+1;
    for i=1:Num_col-1
        spring_num=spring_num+1;
        spring(spring_num).from=part_num;
        spring(spring_num).to=part_num+1;
        spring(spring_num).length=L;
        spring(spring_num).KS=Ks;
        spring(spring_num).KD=Kd;
        part_num=part_num+1;
        Line1(spring_num)=line([x(spring(spring_num).from,1),...
            x(spring(spring_num).to,1)],[x(spring(spring_num).from,2),...
            x(spring(spring_num).to,2)],'marker','o');
%         hold on
    end
end

% spring(28), last horisontal spring

part_num2=0;

for j=1:Num_col
    part_num2=part_num2+1;
    part_num=part_num2;
    for i=1:Num_row-1
        spring_num=spring_num+1;
        spring(spring_num).from=part_num;
        spring(spring_num).to=part_num+Num_col;
        spring(spring_num).length=L;
        spring(spring_num).KS=Ks;
        spring(spring_num).KD=Kd;
        part_num=part_num+Num_col;
        Line1(spring_num)=line([x(spring(spring_num).from,1),...
            x(spring(spring_num).to,1)],[x(spring(spring_num).from,2),...
            x(spring(spring_num).to,2)],'marker','o');
%         hold on
    end
end

% spring(52), last vertical spring

part_num2=0;

for j=1:Num_col-1
    part_num2=part_num2+1;
    part_num=part_num2;
    for i=1:Num_row-1
        spring_num=spring_num+1;
        spring(spring_num).from=part_num;
        spring(spring_num).to=part_num+Num_col+1;
        spring(spring_num).length=L2;
        spring(spring_num).KS=Ks;
        spring(spring_num).KD=Kd;
        part_num=part_num+Num_col;
        Line1(spring_num)=line([x(spring(spring_num).from,1),...
            x(spring(spring_num).to,1)],[x(spring(spring_num).from,2),...
            x(spring(spring_num).to,2)],'marker','o');
%         hold on
    end
end

% spring(73), last right-going diag. spring

part_num2=1;

for j=1:Num_col-1
    part_num2=part_num2+1;
    part_num=part_num2;
    for i=1:Num_row-1
        spring_num=spring_num+1;
        spring(spring_num).from=part_num;
        spring(spring_num).to=part_num+Num_col-1;
        spring(spring_num).length=L2;
        spring(spring_num).KS=Ks;
        spring(spring_num).KD=Kd;
        part_num=part_num+Num_col;
        Line1(spring_num)=line([x(spring(spring_num).from,1),...
            x(spring(spring_num).to,1)],[x(spring(spring_num).from,2),...
            x(spring(spring_num).to,2)],'marker','o');
%         hold on
    end
end

%  spring(94), last left-going diag. spring

for j=1:spring_num
    
        r=x_new(spring(j).from,:)-x_new(spring(j).to,:);
        rdot=v_new(spring(j).from,:)-v_new(spring(j).to,:);
        force=-(spring(j).KS*(norm(r)-spring(j).length)+...
                spring(j).KD*((sum(conj(rdot).*r))/norm(r)))...
                *(r/norm(r));
    
        f(spring(j).from,:)=f(spring(j).from,:)+force;
        f(spring(j).to,:)=f(spring(j).to,:)-force;
end

v_old=v_new-f.*dt/2;

% set_line = @(s) 



set(Line1,'erasemode','xor');

daspect([1,1,1]);
axis([-5,50,-5,10]);

hold on

% my=zeros(Nts,3);

for k=1:Nts
f=zeros(NP,dim);
    for p=1:NP
        for l=1:NumBalls
            dr=(x_old(p,:)-[Pos_ball(l,1),Pos_ball(l,2),0]);
            if norm(dr)<R_ball
            n_hat=dr/norm(dr);
            v_para=dot(v_old(p,:),n_hat)*n_hat;
                v_old(p,:)=v_old(p,:)-2*v_para;
                x_old(p,2)=x_old(p,2)+2*(R_ball-dr(1,2));
            end
        end
    end
   
    
    for j=1:spring_num
        
        r=x_old(spring(j).from,:)-x_old(spring(j).to,:);
        rdot=v_old(spring(j).from,:)-v_old(spring(j).to,:);
        force=-(spring(j).KS*(norm(r)-spring(j).length)+...
                spring(j).KD*((sum(conj(rdot).*r))/norm(r)))...
                *(r/norm(r));
    
        f(spring(j).from,:)=f(spring(j).from,:)+force;
        f(spring(j).to,:)=f(spring(j).to,:)-force;
        fjf(j)=norm(r)-spring(j).length;
    end
    
    f=f+m*g;
    
    v_new=v_old+m\f*dt;
    x_new=x_old+v_new*dt;
    v_k=(v_old+v_new)/2; %mittenvärde
    
    v_old=v_new;
    x_old=x_new;
    
    
%     if k==n*100
     
%       n=n+1;
        Ek(k)=sum(sum((1/2)*(m*(v_k).^2)));
        Ep(k)=sum((-m*g)'*x_new(:,2));
        Efj(k)=sum((1/2)*Ks*fjf.^2);
        Etot(k)=Ek(k)+Ep(k)+Efj(k);
%     end

    v_centmassx(k)=(sum(m*v_new(:,1)))/sum(sum(m));
    v_centmass=(sum(m*v_new))/sum(sum(m));
    
    my(k)=norm((v0-v_centmass)/(gravity*k*dt));

%     for i=1:spring_num
%     
%         set(Line1(i),'xdata',[x_new([spring(i).from],1); x_new([spring(i).to],1)], ...
%                      'ydata',[x_new([spring(i).from],2); x_new([spring(i).to],2)]);
%     end
%      pause(10^(-5));
end

my_av=sum(my)/Nts;

figure
plot(Ek);
hold on
plot(Ep);
hold on
plot(Efj);
hold on
plot(Etot);
legend('Ek','Ep','Efj','Etot')
figure
plot(v_centmassx)
figure
plot(my);

