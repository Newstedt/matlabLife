function [R,dR,dRn]=rotmat2d(th)
%MATMUL2 2D matrix multiplication with Jacobians.
%
%   Y=MATMUL2(R,P) computes the product Y=R*P.
%   [Y,dY]=MATMUL2(R,P) also returns a struct dY with fields dR and dP
%   that contain the analytical Jacobians of Y with respect to R and
%   P, respectively.
%
%   TESTING FUNCTIONS:
%  
%   [Y,dY,dYn]=MATMUL2(R,P) returns a similar struct dYn with the
%   numerical approximations of the Jacobians in dY.
%
%   Y=MATMUL2('selftest') or MATMUL2 selftest runs a small selftest
%   and returns a two-vector with the maximum absolute and relative
%   error between the Jacobians in dY and dYn, as measured by the
%   Frobenius norm.

% Author: Niclas Börlin, niclas.borlin@cs.umu.se
%   2018-11-12: First version.

% Handle selftest.
if ischar(th) && strcmp(th,'selftest'), [R,dR,dRn]=selftest; return; end

% Compute the function.
R = [cos(th), -sin(th); sin(th), cos(th)];

if nargout>1
    % Preallocate the struct for the analytical Jacobians.
    dR=struct('dth',[]);
end

if nargout>2
    % Preallocate the struct for the numerical Jacobians.
    dRn=struct('dth',[]);
end

if nargout>2
    % Compute numerical Jacobians.
    
    % Note: In the definition below, the variable R binds to the parameter
    % used in the call to f. The variable P binds to the variable P
    % known when the function handle was created.
    f=@(th)feval(mfilename,th);
    dRn.dth=jacapprox(f,th);

end

if nargout>1
    % Compute analytical Jacobians.

    % Your code goes here.
    dR.dth = -[sin(th),-cos(th),cos(th),sin(th)]';
    
end

function [v,dR,dRn]=selftest
% Run selftest of the analytical Jacobian. Returns the vector r
% with the absolute and relative errors.

% Set up a reasonable sized problem. Use random numbers to avoid
% selecting a special case. The setup code will be different for each
% residual function.
th=rand(1);

% The test code below can remain the same, except for the parameter list.

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

[~,dR,dRn]=feval(mfilename,th);

fieldNames=fieldnames(dR);
a=-inf;
r=-inf;
for i=1:length(fieldNames)
    fName=fieldNames{i};
    a=max(a,abserr(dR.(fName),dRn.(fName)));
    r=max(r,relerr(dR.(fName),dRn.(fName)));
end

v=[a;r];
