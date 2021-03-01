x=rand(2,1);
A=rand(1e6,2);
[r,J, JJ]=linres(x,A,A*x);

diff = norm(J-JJ,'fro');