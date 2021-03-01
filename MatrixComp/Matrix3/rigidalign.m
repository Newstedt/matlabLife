function [R,t] = rigidalign(P,Q)
%RIGIDALINE - Takes two bodies P and Q. Finds rotation-matrix and
%             translation vector t, for optimal alignment of P with
%             respect to Q.
%
%   MINIMAL WORKING EXAMPLE: We have two 3-point bodies with coordinates
%                            defined by:
%
%           P = [p_x1 p_x2 p_x3; p_y1 p_y2 p_y3]
%           Q = [q_x1 q_x2 q_x3; q_y1 q_y2 q_y3].
%           
%           We want to find 2-by-2 matrix R and 2-by-1 matrix t, that
%           optimally align P with respect to Q.
%
%   >> P = [p_x1 p_x2 p_x3; p_y1 p_y2 p_y3]; %define P
%   >> Q = [q_x1 q_x2 q_x3; q_y1 q_y2 q_y3] %define Q
%   >> [R,t] = rigidalign(P,Q) %find R and t by calling rigidalign

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-12: Initial version .
%
% Function code starts here...

p_cent = [mean(P(1,:));mean(P(2,:))]; %Compute the center of mass fo body P
q_cent = [mean(Q(1,:));mean(Q(2,:))]; %Compute the center of mass fo body Q

X = P-p_cent; %Introduce new normalized variable X for body P
Y = Q-q_cent; %Introduce new normalized variable X for body Q

S = X*Y'; %Set up covariance matrix
[U, sigma, V] = svd(S); %Perform SVD on the convariance matrix
mid = eye(size(V,1)); %Start by def I-matrix to ensure pure rot matrix R
mid(end,end) = det(V*U');%Set last elem on diag of I-matrix=det(V*U')->pure
R = V*mid*U'; %Define R

t = q_cent - R*p_cent; %Set t