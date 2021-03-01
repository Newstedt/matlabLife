function Y=normhomo(X,strip)
%NORMHOMO Normalize homogenous points.
%
%Y=normhomo(X[,strip]);
%X     - homogenous points of any dimension as columns.
%strip - if supplied and non-zero, the trailing one is stripped after
%        normalization.
%Y     - normalized homogenous points (strip==0) or euclidian points
%        (strip~=0).

% v1.0  2002-03-19. Niclas Borlin, niclas@cs.umu.se.

if (nargin<2), strip=0; end

if (strip)
	Y=X(1:end-1,:)./repmat(X(end,:),size(X,1)-1,1);
else
	Y=X./repmat(X(end,:),size(X,1),1);
end
