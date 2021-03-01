function [r,J] = Residual_Function(x,p)
h = Model_Function(x);
r = h-p;
r = reshape(r,[],1);
J = zeros(2*length(x(6:end)),length(x(6:end))+ 5);

for i = 1:length(x(6:end))
    J(2*i-1:2*i,:) = [1 0 cos(x(5))*cos(x(5+i)) -sin(x(5))*sin(x(5+i)) -x(3)*sin(x(5))*cos(x(5+i))-x(4)*cos(x(5))*sin(x(5+i)) zeros(1,length(x(6:end)));...
                        0 1 sin(x(5))*cos(x(5+i)) cos(x(5))*sin(x(5+i)) x(3)*cos(x(5))*cos(x(5+i))-x(4)*sin(x(5))*sin(x(5+i)) zeros(1,length(x(6:end)))];
    J(2*i-1:2*i,i+5) = [-x(3)*cos(x(5))*sin(x(5+i))-x(4)*sin(x(5))*cos(x(5+i));...
        -x(3)*sin(x(5))*sin(x(5+i))+x(4)*cos(x(5))*cos(x(5+i))];
end
end