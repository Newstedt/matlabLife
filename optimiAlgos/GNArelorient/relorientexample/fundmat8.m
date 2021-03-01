function F=fundmat8(X1,X2,normalize,constrain)
%FUNDMAT8 Calculate the fundamental matrix from 8 or more correspondences.
%
%F=fundmat8(X1,X2[,normalize[,constrain])
%X1 - homogenous image coordinates in the first image.
%X2 - homogenous image coordinates in the second image.
%F  - the fundamental matrix such that X2(:,i)'*F*X1(:,i) approx = 0.

if (nargin<3), normalize=true; end
if (nargin<4), constrain=true; end

% Normalize coordinates.
X1=normhomo(X1);
X2=normhomo(X2);

% Calculate normalized transformations.
T1=eye(3);
T2=eye(3);

if (normalize)
    % Translate center of gravity to the origin.
	T1(1:2,3)=-mean(X1(1:2,:),2);
	X1T=T1*X1;
    % Scale such that average point distance to origin is sqrt(2).
	TS1=diag([[1,1]*1./sqrt(mean2(X1T(1:2,:).^2)),1]);
	T1=TS1*T1;

    % Translate center of gravity to the origin.
	T2(1:2,3)=-mean(X2(1:2,:),2);
	X2T=T2*X2;
    % Scale such that average point distance to origin is sqrt(2).
	TS2=diag([[1,1]*1./sqrt(mean2(X2T(1:2,:).^2)),1]);
	T2=TS2*T2;

	% Transform points.
	X1=T1*X1;
	X2=T2*X2;
end

% Extract 2d coordinates.
W1=X1(1:2,:)';
W2=X2(1:2,:)';

% Construct the A matrix.
A=[W2(:,1).*W1(:,1), W2(:,1).*W1(:,2), W2(:,1), ...
   W2(:,2).*W1(:,1), W2(:,2).*W1(:,2), W2(:,2), ...
   W1, ones(size(W1,1),1)];

% Solve for f.
[U,S,V]=svd(A,0);
f=V(:,end);

% Pack into a 3x3 matrix.
F=reshape(f,3,3)';

if (constrain)
    % Find a rank-2 approximation of F.
	[U,S,V]=svd(F);
	S(3,3)=0;
	F=U*S*V';
end

% Calculate fundamental matrix for untransformed points.
F=T2'*F*T1;
