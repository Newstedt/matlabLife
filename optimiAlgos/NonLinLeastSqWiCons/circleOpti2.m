%%
k = 1:1:100;
p = {};

for i = 1:length(k)
    p{i} = circdata(1,9,k(i));
end

figure
hold on
for i = 1:length(k)
    plot(p{i}(1,:),p{i}(2,:),'*')
end
axis equal
%%

[p,C,xNoNoise] = circdata(1,9,0);
x0 = init_x0(p);
x0(4:end) = x0(4:end)+1;

maxIter=50; tol=1e-8;c1=0.1;aMin=1e-3;

fun = @(x) circle_r(x,p);

[x1,code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});

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

[r,J] = circle_r(x1,p);
sigma0 = sqrt(r'*r/(length(p)-3));
c_xx = sigma0^2*inv(J'*J);

sigma_par = diag(c_xx);
sigma_par = sigma_par(1:3);

cov_val = c_xx(1:3,1:3);
corr_val = corrcov(cov_val);

figure
plotcovellipse(x1(1:2),cov_val(1:2,1:2),sqrt(5.99))
axis equal

%%
cK = zeros(length(k),2);
figure
hold on
for i = 1:length(k)
    [p,C,xNoNoise] = circdata(1,9,k(i));
    fun = @(x) circle_r(x,p);
    [x,~,~,~,~] = gaussn(fun,xNoNoise,tol,maxIter,c1,aMin,{});
    cK(i,:) = x(1:2)';
    dtheta = linspace(0,2*pi);
    h = circle_g(x(1:2),x(3),dtheta');
    plot(h(1,:),h(2,:),'-','linewidth', 1.5)
end
axis equal
cK = cK';
cK_mean = mean(cK,2);
cK_cov = cov(cK');

figure
hold on
for i = 1:length(cK)
    scatter(cK(1,i),cK(2,i), 10, 'filled')
end
figure
hold on
plotcovellipse(cK_mean,cK_cov,sqrt(5.99))
plotcovellipse(x1(1:2),cov_val(1:2,1:2),sqrt(5.99))




