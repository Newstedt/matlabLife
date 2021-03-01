function[pos1, pos2] = traj(p1, v1, p2, v2, print)
t = [0:1/60:1];

x1 = p1(1) + v1(1)*t;
x1 = x1';
y1 = p1(2) + v1(2)*t;
y1 = y1';
pos1 = [x1 y1];

x2 = p2(1) + v2(1)*t;
x2 = x2';
y2 = p2(2) + v2(2)*t;
y2 = y2';
pos2 = [x2 y2];

if print == 1
    figure
    hold on
    plot(pos1(:,1), pos1(:,2),'-.','LineWidth',2);
    plot(pos2(:,1), pos2(:,2),'-.','LineWidth',2);
    title('Air crafts trajectories that never cross each other');
end

end