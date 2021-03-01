p = circdata(0,11);
x0 = init_x0(p);
x0(4:end) = x0(4:end)+1;

maxIter=50; tol=1e-8;c1=0.1;aMin=1e-3;

fun = @(x) circle_r(x,p);

[x,code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});

%------------------- Redundancy number ------------------------------------
[r, J, ~] = circle_r(x,p);
U = J*((J'*J)\J');
R = eye(2*length(p)) - U;
RedNum_exp = sum(reshape(diag(R),2,[]),1);
%--------------------------------------------------------------------------

figure
hold on
plot(p(1,:),p(2,:),'*','linewidth' , 1.5)

for i = 1:size(X,2)
    dtheta = linspace(0,2*pi);
    h = circle_g(X(1:2,i),X(3,i),dtheta');
    plot(h(1,:),h(2,:),'-','linewidth', 1.5)
    title([sprintf("\n Gauss-Newton, Iter = %d, code = %d \n", n,code)])
end
legend('Points','1st','2nd','3d')
axis equal
figure
plot(1:length(alphas),alphas)
title("Iteration Number vs Alphas")
xlabel("Iterations")
ylabel("\alpha")

%%
epsR = 1e-3; epsC = 1e-8; nu0 = 0.1; mu = 0.1;
x02 = init_x0_const(p);
con = @(x) circle_gx_const(x,p);
fun2 = @(x) circle_r_const(x,p);
[xc,n,code,l,X2,a,C,L,nus,r,Jc,A]=sqpsq(fun2,con,x02,epsR,epsC,maxIter,{},nu0,mu);

figure
hold on
plot(p(1,:),p(2,:),'*','linewidth' , 1.5)

for i = 1:size(X2,2)
    dtheta = linspace(0,2*pi);
    h2 = circle_g(X2(1:2,i),X2(3,i),dtheta');
    plot(h2(1,:),h2(2,:),'-','linewidth', 1.5)
    title([sprintf("\n SQPSQ, Iter = %d, code = %d \n", n,code)])
end
legend('Points','1st','2nd','3d')
axis equal

%------------------- Redundancy number ------------------------------------
M = [Jc'*Jc, A'; A, zeros(size(A,1))];
N = inv(M);
E = [eye(size(A,2)), zeros(size(A,2),size(A,1))];

P_imp = E*N*E'*Jc';
U = Jc*P_imp;
R_imp = eye(length(U))-U;
sum(reshape(diag(R_imp),2,[]),1)'
% 
% U = J*((J'*J)\J');
% R = eye(2*length(p)) - U;
% RedNum_imp = sum(reshape(diag(R),2,[]),1);
%--------------------------------------------------------------------------
