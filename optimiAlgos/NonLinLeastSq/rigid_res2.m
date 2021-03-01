function [r,J,JJ] = rigid_res2(x,P,Q)
%RIGID_RES2 Compute residual and Jacobians for the lrigid body problem,
%           using chain rule when computing the Jacobian.
%   
%   computes the residual r = vec(R*P+kron(ones(1,length(P)),d)-Q)
%   r = rigid_res2(x,P,Q)
%
%   also return the analytical Jacobian J of the residual:
%   [r,J] = rigid_res2(x,P,Q)
%
%   Minimum working example:
%
%   P = rand(2,5);
%   Q = rand(2,5);
%   x = rand(3,1);
%   [r,J] = rigid_res2(x,P,Q);
%
%   Testing functions:
%
%   further, return a numerical approximation, JJ, of the Jacobian:
%   [r,J,JJ] = rigid_res2(x,P,Q)
%
%   rigid_res2('selftest') runs a small selftest and
%   returns a two-vector with the absolute and relative error between
%   J and JJ, as measured by the Frobenius norm.

% Author: Gustav Nystedt  & Fredrik Gunnarsson,
% guny0007@student.umu.se & frgu0031@student.umu.se
%   2018-11-18: First version.

% Handle selftest.
if ischar(x) && strcmp(x,'selftest'), [r,J,JJ]=selftest; return; end

% Define th, d and R from x
th = x(1);
d = x(2:end);
R = [cos(x(1)), -sin(x(1)); sin(x(1)), cos(x(1))];

% Compute the residual.
r = vec(R*P+kron(ones(1,length(P)),d)-Q);

% Get dR and dY = d(PQ) for further use in chain rule
[~,dR,~] = rotmat2d(th);
[~,dY,~] = matmul2(R,P);

if nargout>2
    % Return numerical Jacobian.
    f=@(x)feval(mfilename,x,P,Q);
    JJ=jacapprox(f,x);
end

if nargout>1
    % Return analytical Jacobian, using chain rule in first column
    J = [dY.dR*dR.dth, kron(ones(length(P),1),eye(2))];
end
end

function [r,J,JJ]=selftest
% Run selftest of the analytical Jacobian. Returns the vector r
% with the absolute and relative errors.

% Set up a reasonable sized problem. Use random numbers to avoid
% selecting a special case. The setup code will be different for each
% residual function.
m=2;
n=5;
P=rand(m,n);
Q=rand(m,n);
x=rand(3,1);

% The test code below can remain the same, except for the parameter list.

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

[~,J,JJ]=feval(mfilename,x,P,Q);

r=[abserr(J,JJ);relerr(J,JJ)];
end

function v=vec(x)
% Specify vec here to make jacapprox self-contained.

v=x(:);
end

