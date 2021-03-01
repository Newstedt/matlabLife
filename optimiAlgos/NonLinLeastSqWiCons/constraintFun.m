function [c, ceq] = constraintFun(theta)

dim = floor(sqrt(length(theta)));
for i = 1:dim
    const(i) = sum(theta( (i-1)*dim+3 : i*dim+2 ).^2) - 1;
end

c = [];
ceq = theta + sum(const);
