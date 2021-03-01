% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D diffusion problem for general boundary conditions
% This code can be used for a flux or set temperature boundary condition
% depending on a flag. It is shown for demonstration how to make your codes
% more sophisticated.
% This is a line by line model solver which is somewhat advanced
% Here we solve
% u^(n+1)=u^n+dt(a*u^n+source+boundary) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for N=16

%N=32;                        % Number of volumes
Uinitial=0;               % initial temperature distribution (T=U)
v=0.000217;                       % Diffusivity (k=v)
h=0.3e-2;                     % Length (L=h)
dy=h/N;  
dt=0.25*dy^2/v;                       % Time step
dy=h/N;                     % lenght of cell
dyConst=h/32;
time_final=0.1;             % Final time
d=v/dyConst^2;            % Define new function
Ub = 0;                     % Value of T on boundary
U0 = [0 10];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate in time system of equations using Forward Euler
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uold=Uinitial*ones(N,1);
timeCount=0;
stagCount=0;
for time=dt:dt:time_final
% Update boundary conditions. Must be done here for time varying boundary 
% conditions
    timeCount=timeCount+1;
    i=1;
    if(time == dt)
        U0_use=U0(1);
    
    else
        U0_use=U0(2);
    end
        
    unew(i)=uold(i)+dt*d*(uold(i+1)-3*uold(i)+2*U0_use);
        
    for i=2:1:N-1
        unew(i)=uold(i)+dt*d*(uold(i+1)-2*uold(i)+uold(i-1));
    end
    
    i=N;
    unew(i)=uold(i)+dt*d*(2*Ub-3*uold(i)+uold(i-1));

    % Exchange old and new
    uold=unew;
    uHalf(timeCount)=unew(floor(length(unew)/2));
%         if(time>100*dt)
%             if(uHalf(timeCount)/uHalf(timeCount-100)<1.05)
%                 stagCount=stagCount+1;
%                 stagTime(stagCount)=time;
%             end
%         end
    %ysol=dy/2:dy:h-dy/2;
    %ymid(timeCount)=ysol(floor(length(ysol)/2));
    timePlot(timeCount)=time;
    
end
plot(timePlot,uHalf)
hold on
end
%figure;plot(timePlot,uHalf)
% ysol=dy/2:dy:h-dy/2;
% stagTime(1)
% figure;plot(unew,ysol,'-bs') 
