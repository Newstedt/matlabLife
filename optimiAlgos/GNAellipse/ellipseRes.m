function [r,J] = ellipseRes(x,p)

h = ellipseMod(x);
nabla1 = zeros(2*(length(x)-5),5);
nabla2 = zeros(2*(length(x)-5),length(x)-5);
r = h-p;
r=reshape(r,[],1);

for i = 1:length(x)-5
   
   nabla1(2*i-1:2*i,:) = ...
   [1 0 cos(x(5))*cos(x(5+i)) -sin(x(5))*sin(x(5+i)) -x(3)*sin(x(5))*cos(x(5+i))-x(4)*cos(x(5))*sin(x(5+i));
    0 1 sin(x(5))*cos(x(5+i)) cos(x(5))*sin(x(5+i)) x(3)*cos(x(5))*cos(x(5+i))-x(4)*sin(x(5))*sin(x(5+i))];
    
   nabla2(2*i-1:2*i,i) = ...
   [-x(3)*cos(x(5))*sin(x(5+i))-x(4)*sin(x(5))*cos(x(5+i));
    -x(3)*sin(x(5))*sin(x(5+i))+x(4)*cos(x(5))*cos(x(5+i))];
end

nabla = [nabla1,nabla2];
J = nabla;



