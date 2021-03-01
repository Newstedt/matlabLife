function Y=euclidian(X)
%EUCLIDIAN Convert homogenous coordinates to Euclidian.
%
%Y=euclidian(X)
%X - (k+1)-by-n matrix of homogenous k-dimensional coordinates.
%Y - k-by-n matrix of Euclidian k-dimensional coordinates.

Y=X(1:end-1,:)./X(repmat(end,1,size(X,1)-1),:);
