function Y=homogenous(X)
%HOMOGENOUS Convert Euclidian coordinates to homogenous.
%
%Y=homogenous(X)
%X - k-by-n matrix of Euclidian k-dimensional coordinates.
%Y - (k+1)-by-n matrix of homogenous k-dimensional coordinates.

Y=[X;ones(1,size(X,2))];
