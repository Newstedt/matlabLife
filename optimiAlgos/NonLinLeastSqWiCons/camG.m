function [c,A,AA] = camG(x,p)
%camG - Function that setup constraint function  and its 
%       Jacobian for photogrammetry problem.
%
%   Minimum working example:
%
%   Used inside sqpsq.m:
%   [c,A,AA] = camG(x,p)

% Author: Gustav Nystedt, guny0007@ad.umu.se
%         Fredrik Gunnarsson, frgu0031@ad.umu.se
% 2018-12-20: First version.

%Unpack x
R = x(1:27);
d = x(28:36);
q = reshape(x(37:end),3,[]);

% split R
R1 = reshape(R(1:9),[],3)';
R2 = reshape(R(10:18),[],3)';
R3 = reshape(R(19:end),[],3)';

% redefine R
R = {R1,R2,R3};

% define components of constraint function
g1 = reshape((R1-eye(3))',[],1);
g2 = d(1:3);
g3 = [R2(:,1)'*R2(:,2); R2(:,1)'*R2(:,3); R2(:,2)'*R2(:,3)];
g4 = [R2(:,1)'*R2(:,1)-1; R2(:,2)'*R2(:,2)-1; R2(:,3)'*R2(:,3)-1];
g5 = [R3(:,1)'*R3(:,2); R3(:,1)'*R3(:,3); R3(:,2)'*R3(:,3)];
g6 = [R3(:,1)'*R3(:,1)-1; R3(:,2)'*R3(:,2)-1; R3(:,3)'*R3(:,3)-1];

% set constraint function
c = [g1;g2;g3;g4;g5;g6];

if nargout>2 %We want the numerical jacobian
    f = @(x) camG(x,p);
    AA = jacapprox(f,x);
end

if nargout>1 %We want the analytical jacobian
   A_g1 = eye(9);
   A_g2 = eye(3);
   j = [2,3,1];
   
   for k = 1:2
       col = 1;
       for i = 1:9
           A_g3Temp(:,i) = [R{k+1}(col,j(1)), R{k+1}(col,j(2)), R{k+1}(col,j(3))]';
           j = [j(end),j(1:2)];
           if mod(i,3) == 0
               col = col+1;
           end
       end
       era = ones(3)-fliplr(eye(3));
       A_g3{k} = A_g3Temp.*[era,era,era];
   end
   
   for k = 1:2
       for i = 1:3
           A_g4{k}(1:3,3*(i-1)+1:3*i) = 2*diag(R{k+1}(i,:));
       end
   end
   A = zeros(24,57);
   A(1:9,1:9) = A_g1;
   A(10:12,28:30) = A_g2;
   A(13:18,10:18) = [A_g3{1};A_g4{1}];
   A(19:24,19:27) = [A_g3{2};A_g4{2}];
end

