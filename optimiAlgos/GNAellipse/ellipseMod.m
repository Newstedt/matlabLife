function [h] = ellipseMod(x)

h = zeros(length(x)-5,2);

for i = 1:length(x)-5
    h(i,:) = x(1:2) + [cos(x(5)), -sin(x(5)); sin(x(5)), cos(x(5))]*...
            (x(3:4).*eye(2))*[cos(x(i+5));sin(x(i+5))];
end
h = h';