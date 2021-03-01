function plot_func(X,P,alphas,n,code,title_)
% plot_func - Function that plots the disired optimization plots  
%
%   Minimum working example: Used inside Test
%                           >> plot_func(X,P,alphas,n,code,title_)
%                               
% 
%   X = Data from optimization 
%   P = Given points
%   alphas = alphas from optimization algorithm
%   n = number of iteration computed
%   code = sucess or failure from optimization
%   title_ = title for the plots

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

plot(P(1,:),P(2,:),'*','linewidth' , 1.5)
hold all

for i = 1:size(X,2)
    dtheta = linspace(0,2*pi);
    h = circle_g(X(1:2,i),X(3,i),dtheta');
    plot(h(1,:),h(2,:),'-','linewidth', 1.5)
    title([sprintf("\n Iter = %d, code = %d \n %s",n,code,title_)])
end
legend('Points','1st','2nd')
axis equal
if alphas ~= 0
    figure
    plot(1:length(alphas),alphas)
    title("Iteration Number vs Alphas")
    xlabel("Iterations")
    ylabel("\alpha")
end

%% residual
if alphas==0
    figure
    for i = 1:size(X,2)
        [r,J] = circle_r_const(X(:,i),P);
        res(i)= r'*r/2;
    end
    semilogy(0:length(res)-1,res)
    xlabel("Iteration")
    ylabel("Objective Function")
end

