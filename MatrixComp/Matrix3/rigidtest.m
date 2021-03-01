% Load the house coordinates.
P=house;
% Pick an angle between 0..90 degrees
alpha=rand*pi/2
% Create a rotation matrix
R0=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
% Verify that R0 is a rotation matrix.
det(R0) % (should be close to unity)
% Generate a random shift.
m=10+rand(2,1);
% Rotate and shift the house.
PT=R0*P+repmat(m,1,size(P,2));
% Perturb the transformed house with noise.
noiseLevel=0.25;
Q=PT+randn(size(PT))*noiseLevel;
% Plot the original and perturbed house.
plot(P(1,:),P(2,:),'bo-',Q(1,:),Q(2,:),'gx-'); axis equal
% At this point, only P and Q are assumed to be know.
% Now, compute the rigid-body transformation between P and Q...
[R,t]=rigidalign(P,Q);
% Apply the inverse operation to Q to align with P.
PA=R'*(Q-repmat(t,1,size(Q,2)));
% Plot the re-aligned house.
line(PA(1,:),PA(2,:),'marker','*','linestyle','-','color','r');