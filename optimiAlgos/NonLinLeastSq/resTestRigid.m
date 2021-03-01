%RESTESTRIGID Script that runs tests for answering question 1.5
%
%   Minimum working example:
%
%   >> resTestRigid

% Author: Gustav Nystedt          & Fredrik Gunnarsson,
%         guny0007@student.umu.se & frgu0031@student.umu.se
%   2018-11-18: First version.

% Define P and Q from exercise
P = [1 3 3;1 1 2];
Q = [-1 1 1;2 3 4];

% Visualise P and Q
figure
hold on
scatter(Q(1,:),Q(2,:))
scatter(P(1,:),P(2,:))

% Define x as given in exercise
x1 = [0.3488, -1.4039, 0.9495]';

% Use rigid_res2 to get jacobian and residuals for x, P and Q
[r1,J1,JJ1] = rigid_res2(x1,P,Q);

% Define variable range for testing surrounding points
theta = 0.3478:0.0001:0.3498;
d1 = -1.4050:0.0001:-1.403;
d2 = 0.9485:0.0001:0.9505;

l = 1; % Define stepping variable
nr = zeros(length(theta)*length(d1)*length(d2),1); % Pre-allocation
x = cell(length(theta)*length(d1)*length(d2),1); % Pre-allocation

% Loop over variables to test surrounding points
for i = 1:length(theta)
    for j = 1:length(d1)
        for k = 1:length(d2)
            x{l} = [theta(i), d1(j), d2(k)]'; % Extend x
            [r,~,~] = rigid_res2(x{l},P,Q); % Call rigid_res2 to get res
            nr(l) = norm(r); % Compute norm of residual and save to vector
            l = l + 1; % Increase stepping variable
        end
    end
end

% Plot residual norms
figure
plot(nr)

% Test to check if rigidalign produces same result as in previouse exercise
[R,t]=rigidalign(P,Q);
xApp = [acos(R(1,1));t(1);t(2)];

