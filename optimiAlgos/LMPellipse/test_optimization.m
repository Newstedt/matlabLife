function [x,code,n,X,dk,rhos,stepType] = test_optimization(Image,Type,Initial_Type)
figure
maxIter=50; tol=1e-8;mu=0.25;eta=0.75;c1=0.1;aMin=1e-3;

% Get datapoints
p = ellipse_data(Image(1),Image(2));

% Define function
fun = @(y) Residual_Function(y,p);

% compute X0
if (strcmp(Initial_Type,'ellip_x0'))
    [c,ax,delta,theta] = ellip_x0(p(1,:)',p(2,:)');
else
    [c,ax,delta,theta] = ellip_x0_simple(p(1,:)',p(2,:)');

end
x0 = [c',ax',delta',theta]';
h = Model_Function(x0);
%%
if strcmp(Type,'Gauss')
    [x,code,n,X,alphas] = gaussn(fun,x0,tol,maxIter,c1,aMin,{});
    hold all
    plot(p(1,:),p(2,:),'*','linewidth' , 1.5)
    for i = 1:size(X,2)
        dtheta=linspace(0,2*pi);
        h = Model_Function([X(1:5,i);dtheta']);
        plot(h(1,:),h(2,:),'-','linewidth', 1.5)
        title([sprintf("\n Gauss-Newton, Iter = %d, code = %d, image: [%d %d] \n", n,code,Image(1),Image(2))])

    end
    legend('Points','1st','2nd','3d')
axis equal
end

%% Test Levenberg
if strcmp(Type,'LM')
    d0 = norm(x0);
    [x,code,n,X,dk,rhos,stepType]=levenberg(fun,x0,tol,maxIter,d0,mu,eta,{});
    hold all
    plot(p(1,:),p(2,:),'*','linewidth' , 1.5)
    for i = 1:size(X,2)
        dtheta=linspace(0,2*pi);
        h = Model_Function([X(1:5,i);dtheta']);
        plot(h(1,:),h(2,:),'-','linewidth', 1.5)
        title([sprintf("\n Levenberg, Iter = %d, code = %d, image: [%d %d] \n", n,code,Image(1),Image(2))])
    end
    axis equal
    legend('Points','1st','2nd','3d')
end


%% residual
figure
for i = 1:size(X,2)
    h = Model_Function(X(:,i));
    r = h - p;
    r = reshape(r,[],1);
    res(i)= r'*r/2;
    title([sprintf("\n Levenberg, Iter = %d, code = %d, image: [%d %d] \n", n,code,Image(1),Image(2))])
end
plot(1:length(res),res)
title([sprintf("\n %s, Iter = %d, code = %d, image: [%d %d] \n",Type, n,code,Image(1),Image(2))])

