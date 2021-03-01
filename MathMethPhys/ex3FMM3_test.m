clc
clear all
%Number of points, dimension, number of time steps and delta t
NP=32; dim=3; Nts=1500; dt=0.002; Num_col=8; Num_row=4; r=zeros(1,dim);
rdot=zeros(1,dim); f=zeros(NP,dim); gravity=1; g=zeros(NP,dim); 
g(:,2)=-gravity;

%Properties of string
NS=96; L=1; L2=sqrt(2); Ks=100; Kd=0;

%initial values
mval=ones(32,1); m=diag(mval);

%initial values balls
NumBalls=16;
X_ball=1;
Y_ball=-2;
R_ball=L;

%position-vector first ball
Pos_ball=zeros(NumBalls,dim);
Pos_ball(1,:)=[X_ball, Y_ball, 0];

%placing balls

for i=2:NumBalls
    Pos_ball(i,:)=Pos_ball(i-1,:)+[R_ball+0.1*R_ball*randn,0,0];
    BOLL_a=rectangle('Position',[Pos_ball(i,1)-R_ball,Pos_ball(i,2)-R_ball,2*R_ball,2*R_ball],...
    'Curvature',[1,1],'EdgeColor','r','FaceColor','w');
end

spring_num=0;
k=0;

for i=1:Num_row
    for j=1:Num_col
        k=k+1;
        x(k,:)=[j,i,0];
        v(k,:)=[6,-10,0];
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
axis([-5,25,-5,15]);

hold on

for k=1:Nts
f=zeros(NP,dim);
    for p=1:NP
        for l=1:NumBalls
            t=norm((x_new(p,:)-[Pos_ball(l,1),Pos_ball(l,2),0]));
            if t<R_ball
                v_old(p,:)=v_old(p,:).*[1,0,1];
            end
        end
    end
    
    v_new=v_old;
    
    for j=1:spring_num
        
        r=x_new(spring(j).from,:)-x_new(spring(j).to,:);
        rdot=v_new(spring(j).from,:)-v_new(spring(j).to,:);
        force=-(spring(j).KS*(norm(r)-spring(j).length)+...
                spring(j).KD*((sum(conj(rdot).*r))/norm(r)))...
                *(r/norm(r));
    
        f(spring(j).from,:)=f(spring(j).from,:)+force;
        f(spring(j).to,:)=f(spring(j).to,:)-force;
    end
    
    f=f+m*g;
    
    v_new=v_old+m\f*dt;
    x_new=x_old+v_new*dt;
    v_k=(v_old+v_new)/2; %mittenvärde
    
    v_old=v_new;
    x_old=x_new;
    
    for i=1:spring_num
    
        set(Line1(i),'xdata',[x_new([spring(i).from],1); x_new([spring(i).to],1)], ...
                     'ydata',[x_new([spring(i).from],2); x_new([spring(i).to],2)]);
    end
    pause(0.0001);
end

