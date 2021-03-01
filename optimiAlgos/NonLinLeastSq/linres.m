function [r,J,JJ]=linres(x,A,b)
%LINRES Compute residual and Jacobians for the linear problem.
%
%   r=LINRES(x,A,b) computes the residual r=A*x-b.
%
%   [r,J]=LINRES(x,A,b) also returns the analytical Jacobian J of the residual.
%
%   Minimum working example:
%
%   x=rand(2,1);
%   A=rand(3,2);
%   [r,J]=LINRES(x,A,A*x);
%
%   Testing functions:
%
%   [r,J,JJ]=LINRES(x,A,b) furthermore returns a numerical approximation
%   JJ of the analytical Jacobian.
%
%   LINRES('selftest') or LINRES selftest runs a small selftest and
%   returns a two-vector with the absolute and relative error between
%   J and JJ, as measured by the Frobenius norm.

% Author: Niclas Börlin, niclas.borlin@cs.umu.se
%   2018-11-08: First version.

% Handle selftest.
if ischar(x) && strcmp(x,'selftest'), [r,J,JJ]=selftest; return; end

% Compute the residual.
r = A*x - b; % Your stuff goes here.

if nargout>2
    % Return numerical Jacobian.
    f=@(x)feval(mfilename,x,A,b);
    JJ=jacapprox(f,x);
end

if nargout>1
    % Return analytical Jacobian.
    J = A; % Your stuff goes here.
end


function [r,J,JJ]=selftest
% Run selftest of the analytical Jacobian. Returns the vector r
% with the absolute and relative errors.

% Set up a reasonable sized problem. Use random numbers to avoid
% selecting a special case. The setup code will be different for each
% residual function.
m=4;
n=3;
A=rand(m,n);
b=rand(m,1);
x=rand(n,1);

% The test code below can remain the same, except for the parameter list.

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

[~,J,JJ]=feval(mfilename,x,A,b);

r=[abserr(J,JJ);relerr(J,JJ)];