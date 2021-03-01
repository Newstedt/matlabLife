%%
P = [0 1 2;0 1 0];
Q = [3 4 3;3 2 1];

figure
hold on
scatter(Q(1,:),Q(2,:))
scatter(P(1,:),P(2,:))

[r,J,JJ] = rigid_res2([-pi/2,3,3]',P,Q);

%%
rigid_res1('selftest')
