%%
R = linspace(1,3,27);
d = linspace(1,2,9);
q = linspace(1,8,21);

R1 = reshape(R(1:9),[],3)';
R2 = reshape(R(10:18),[],3)';
R3 = reshape(R(19:end),[],3)';

x = [R,d,q]';
p = randn(3,21);

[r,J,JJ] = camera_r(x,p);
[gx,A,AA] = camG(x,p);
AA = full(AA);
error = AA-A;
