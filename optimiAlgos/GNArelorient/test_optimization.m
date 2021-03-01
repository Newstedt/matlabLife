function [x,code,n,X,dk,rhos,stepType,p] = test_optimization(Image,Type,Initial_Type)
%TEST_OPTIMIZATION Script which compute and plot results from Test.m
%
%   Minimum working example:
%   
%   >> Test

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-05: First version.

figure
maxIter=50; tol=1e-8;mu=0.25;eta=0.75;c1=0.1;aMin=1e-3;dk=0;rhos=0;
stepType=0;

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
figure
plot(1:length(alphas),alphas)
title("Iteration Number vs Alphas")
xlabel("Iterations")
ylabel("\alpha")
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
    
    
    figure
    stem(1:length(rhos),rhos)
    xlabel("Iterations")
    ylabel("\rho_k")
    %ylim([0 1])

    
    figure
    semilogy(1:length(dk),dk)
    xlabel("Iterations")
    ylabel("\Delta_k")
    
    figure
    plot(1:length(stepType),stepType,'*')
    xlabel("Iterations")
    ylabel("Step type")
    ylim([-0.2 2.2])

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
semilogy(0:length(res)-1,res)
title([sprintf("\n %s, Iter = %d, code = %d, image: [%d %d] \n",Type, n,code,Image(1),Image(2))])
xlabel("Iteration")
ylabel("Objective Function")

