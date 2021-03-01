clear all
clc
%%
c1 = 0.1;
aMin = 1e-3;
mu = 0.25;
eta = 0.75;
tol = 1e-8;
maxIter = 50;

p=ellipse_data(2,2);
[c,ax,delta,theta]=ellip_x0(p(1,:)',p(2,:)');

x(1) = c(1);
x(2) = c(2);
x(3) = ax(1);
x(4) = ax(2);
x(5) = delta;

x=x';
d0 = norm(x);

for i = 1:length(theta)
    x(i+5) = theta(i);
end

[r,J] = ellipseRes(x,p);

f = @(x) ellipseRes(x,p);

%%
JJ = jacapprox(f,x);
JJ = full(JJ);

[x2,code,n,X,alphas] = gaussn(f,x,tol,maxIter,c1,aMin,{});
h = {};

% figure
% hold on 
% axis equal
% scatter(p(1,:),p(2,:));
for i = 1:size(X,2)
    h{i} = ellipseMod(X(:,i));
    %scatter(h{i}(1,:),h{i}(2,:))
    plotellipse_with_res(gca,X(:,i),p);
    axis equal
end


