function [P1,P2,X,inFront,sp]=camsfrome(E,x,xp)
%CAMSFROME Calculate cameras from an essential matrix.
%
%[P1,P2,X,inFront,sp]=camsfrome(E,x,xp)
%E       - 3x3 essential matrix.
%x       - 3xn matrix with normalized 2d points in image 1.
%xp      - 3xn matrix with normalized 2d points in image 2.
%P1      - Canonic camera [I,0].
%P2      - Other camera [R,t] with baseline length |t|=1.
%X       - 4xn matrix of reconstructed 3d points with respect to P1.
%inFront - logical n-vector specifying if point X(:,i) is in front of
%          both cameras.
%sp      - sorted vector with number of points in front of each camera config.

[U,S,V]=svd(E);
if (det(U)<0)
	U(:,3)=-U(:,3);
end
if (det(V)<0)
	V(:,3)=-V(:,3);
end

W=[0,-1,0;1,0,0;0,0,1];
Z=[0,1,0;-1,0,0;0,0,0];

% Camera 1.
P1=[eye(3),zeros(3,1)];

% Four choices for camera 2
C2={[U*W*V',U(:,3)],[U*W*V',-U(:,3)],[U*W'*V',U(:,3)],[U*W'*V',-U(:,3)]};
% Reconstructed points for each choice.
XX=cell(size(C2));
% Is the z coordinate positive?
zp=zeros(length(C2),size(x,2));
for j=1:length(C2)
	P2=C2{j};
	% Reconstruct points.
    X=zeros(4,size(x,2));
    for i=1:size(x,2)
        X(:,i)=pm_forwintersect1([P1;P2],[x(:,i),xp(:,i)]);
    end
	% Check what points are in front of both cameras.
	d1=ptdepth(P1,X);
	d2=ptdepth(P2,X);
	zp(j,:)=d1<0 & d2<0;
    XX{j}=X;
end

% Select the camera pair that has the largest number of points in front
% of both cameras.
[sp,ii]=sort(-sum(zp,2));
sp=-sp;
P2=C2{ii(1)};
inFront=zp(ii(1),:);
X=XX{ii(1)};
