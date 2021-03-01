function [Y,dY,dYn]=matmul2(R,P)
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
if ischar(R) && strcmp(R,'selftest'), [Y,dY,dYn]=selftest; return; end

% Compute the function.
Y = R*P; % Your code goes here.

if nargout>1
    % Preallocate the struct for the analytical Jacobians.
    dY=struct('dR',[],'dP',[]);
end

if nargout>2
    % Preallocate the struct for the numerical Jacobians.
    dYn=struct('dR',[],'dP',[]);
end

if nargout>2
    % Compute numerical Jacobians.
    
    % Note: In the definition below, the variable R binds to the parameter
    % used in the call to f. The variable P binds to the variable P
    % known when the function handle was created.
    f=@(R)feval(mfilename,reshape(R,2,2),P);
    dYn.dR=jacapprox(f,R);

    % Note: In the definition below, the variable P binds to the parameter
    % used in the call to f. The variable R binds to the variable R
    % known when the function handle was created.
    f=@(P)feval(mfilename,R,reshape(P,2,[]));
    dYn.dP=jacapprox(f,P);
end

if nargout>1
    % Compute analytical Jacobians.

    % Your code goes here.
    dY.dR = kron(P',eye(2));
    dY.dP = kron(eye(length(P)),R);
end

function [v,dY,dYn]=selftest
% Run selftest of the analytical Jacobian. Returns the vector r
% with the absolute and relative errors.

% Set up a reasonable sized problem. Use random numbers to avoid
% selecting a special case. The setup code will be different for each
% residual function.
n=4;
R=rand(2);
P=rand(2,n);

% The test code below can remain the same, except for the parameter list.

abserr=@(A,B)norm(A-B,'fro');
relerr=@(A,B)abserr(A,B)/norm(A,'fro');

[~,dY,dYn]=feval(mfilename,R,P);

fieldNames=fieldnames(dY);
a=-inf;
r=-inf;
for i=1:length(fieldNames)
    fName=fieldNames{i};
    a=max(a,abserr(dY.(fName),dYn.(fName)));
    r=max(r,relerr(dY.(fName),dYn.(fName)));
end

v=[a;r];
