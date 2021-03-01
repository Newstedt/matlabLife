%%
A = [1 2 3; 0 5 6; 0 0 9]; a = [1.3; 4.2; 12.2];

B = [2 4 1; 0 4 7; 0 0 2]; b = [2.1; 3.4; 1.3];

C = [87 23 54; 0 28 51; 0 0 91]; c = [3.88; 16; 9.5];

D = [1 2 3; 4 5 6];
d = [3 4];

E = [1 1 1; 0 1 1; 0 0 1];
e = [1 1 1];

nrm.a = norm(backsubst(A,a) - A\a);
nrm.b = norm(backsubst(B,b) - B\b);
nrm.c = norm(backsubst(C,c) - C\c);

%%
p1 = [0; 0]; p2 = [0; -1];
v1 = [1; 1]; v2 = [-1; 1];

[x, y, t1, t2] = collision(p1, v1, p2, v2);

%%
p1 = [0; 0]; p2 = [1; 0];
v1 = [1; 1]; v2 = [-1; 1];

pos1 = traj(p1, v1);
pos2 = traj(p2, v2);

figure
hold on
plot(pos1(:,1), pos1(:,2));
plot(pos2(:,1), pos2(:,2));

collision()

function[pos] = traj(p, v)
t = [0:1:10];

x = p(1) + v(1)*t;
x = x';
y = p(2) + v(2)*t;
y = y';
pos = [x y];

end



