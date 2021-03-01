function test_optimization(Image)

maxIter=250;
tol=1e-8;
mu=0.25;
eta=0.75;
c1=0.1;
aMin=1e-3;
p = ellipse_data(Image(1),Image(2));
% Define function
fun = @(y) Residual_Function(y,p);

% compute X0
[c,ax,delta,theta] = ellip_x0(p(1,:)',p(2,:)');

% Perturbation
x0 = [c'+c'*0.1,ax'-ax'.*[0.2,0.7],delta'+delta*10,theta]';

h = Model_Function(x0);
%% check accuracy of jacobian 

[r,J] = Residual_Function(x0,p);
JJ = jacapprox(fun,x0);

Error = max(max(abs(JJ-J)));
fprintf('The largest absolut error for the Jacobian is %.2d \n', Error)