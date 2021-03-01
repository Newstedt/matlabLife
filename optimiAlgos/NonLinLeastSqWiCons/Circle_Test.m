% Circle_Test - Function that run the optimizition for circle parameters 
%
%   Minimum working example:
%                           >> Circle_Test
%

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.


maxIter=50; tol=1e-8; mu = 0.1; c1=0.1;aMin=1e-3;
stepType=0; nu0=0.1; epsC = 1e-8; epsR=1e-3;

P = circdata(0,11); 
x0 = init_x0(P);
x0(4:end) = x0(4:end)+1;

%% No constraints
fun = @(y) circle_r(y,P);
[x,code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});
title="Explicit Parametrization";
plot_func(X,P,alphas,n,code,title)

%% 3.1.4 Redundancy Explicit
m = length(x(4:end));
[r,JJ,J]=circle_r(x,P);
R = eye(2*m) -J*((J'*J)\J'); 
R = sum(reshape(diag(R),2,[]),1)';

%%
fprintf(" Problem & $c_x$ & $c_y$ & $r$ \\\\ \n")
fprintf("%s & %.3d & %.3d & %.3d\\\\ \n",title,x(1),x(2),x(3))

%% Constraint problem
x0 = init_x0_c(P);
fun = @(y) circle_r_const(y,P);
con = @(y) circle_gx_const(y,P);
[x,n,code,l,X,a,C,L,nus,r,J,A]...
    =sqpsq(fun,con,x0,epsR,epsC,maxIter,{},nu0,mu);
title = "Implicit Parametrization";
plot_func(X,P,a,n,code,"Implicit Parametrization");
fprintf("%s & %.3d & %.3d & %.3d\\\\ \n",title,x(1),x(2),x(3))

%% 3.1.4 Redundancy Implicit
E = [eye(length(J)),zeros(size(A,2),size(A,1))];
N = [J'*J A';A zeros(size(A,1))]^-1;
P = E*N*E'*J';
U = J*P;
R1=eye(length(U))-U;
R2 = sum(reshape(diag(R1),2,[]),1)';
ColumnTitle =["p_1","p_2","p_3","p_4","p_5","p_6","p_7","p_7","p_9","p_10","p_11"];
%To_Latex([R,R2],ColumnTitle,["Explicit","Implicit"])

figure
hold all
plot(R2,'-*','linewidth',1.5)
plot(R,'-x','linewidth',1.5)
legend("Implicit Parameterization","Explicit Parameterization")

%% 3.1.5
    figure
    hold all
for k=1:100
    p = circdata(1,9,k);
    plot(p(1,:),p(2,:),'*')
    axis equal
end

%% 3.1.5, investigation 2.1
Parameters = ["\\sigma_0$","c_x","c_y","r"];
[P,cc,xx] = circdata(1,9,0);
x0 = init_x0(P);
x0(4:end) = x0(4:end)+1;
fun = @(y) circle_r(y,P);
[x,code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});
title="Explicit Parametrization";
plot_func(X,P,alphas,n,code,title)
fprintf("\n\n\n\n Problem & $c_x$ & $c_y$ & $r$ \\\\ \n")
fprintf("%s & %.3f & %.3f & %.3f\\\\ \n",title,x(1),x(2),x(3))
[r,J] = circle_r(x,P);
redundancy=2;
sigma(1) = sqrt(r'*r/(length(P)-3));
c_xx = sigma(1)^2*inv(J'*J);

for i = 1:3
    parameter_sigma(1,i) = sqrt(c_xx(i,i));
end

fprintf("\n\n\n& %s & %s & %s & %s \\\\ \n",Parameters)
for i =1:1
        fprintf("[%.2d] & [%.2d] & [%.2d]& [%.2d] \\\\ \n",[sigma(i),parameter_sigma(i,:)])
end
Parameters = ["Image§","\\sigma_0$","c_x","c_y","r"];
RowTitle=["[1,9,0]"];
%To_Latex([sigma,parameter_sigma],RowTitle,Parameters)


%% Plot cov ellipse
figure 
xx=x(1:2);
plotcovellipse(xx,c_xx(1:2,1:2),sqrt(5.99))
hold on
%figure
plot(P(1,:),P(2,:),'*')
legend("Uncertainty Ellips","Circle Center","Points")
axis equal
%% 3.1.5 Simulation
   hold all
   r=0.1063;
   x=[];
for k=1:100
    [p,C,x0]=circdata(1,9,k);
    fun = @(y) circle_r(y,p);
    [x(k,:),code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});
    plot(X(1,end),X(2,end),'*')
end
format long
%%3.1.5-2
Mean = mean(x(:,1:2)',2);
Covariance = cov(x(:,1:2));

%%3.1.5-3
figure
plotcovellipse(xx,c_xx(1:2,1:2),sqrt(5.99));
plotcovellipse(Mean,Covariance,sqrt(5.99));

