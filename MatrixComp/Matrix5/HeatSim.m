%HEATSIM - Script that solves heat equation in 2d-plate, using pcg that
%          uses ichol for preconditioning.
%
%   MINIMAL WORKING EXAMPLE: 
%
%   >> HeatSim;

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-11-01: Initial version .
%
% Function code starts here...

% HeatSim  Similation of heat in 2D plate
tic
% Set number of intervals in spatial direction
N=400; M=N-1;

% Compute spatial stepsize
h=1/N;

% Define grid
t=linspace(1,M,M)*h; [x, y]=meshgrid(t,t);

% Define initial temperatur distribution
u=g(x,y);

% Suppress details from plot to make it readable
nb=floor(N/40); tau=t(1:nb:end,1:nb:end);

% Plot the initial distribution of heat
u1=u(1:nb:end,1:nb:end); surf(tau,tau,u1);

% Maximum temperature display
maxtemp=200;

% Enforce axis
axis([0 1 0 1 0 maxtemp]);

% Hold graphics
hold; 

% Set the simulation time
T=0.03;

% Set the number of timesteps
count=120;

% Compute time step
k=T/count;

% Compute critical parameter lambda
lambda=k/h^2; %multiply by 0.025 to get 1?

% Generate the matrix
A=Heat(N,lambda);

% Generate right hand side
b=reshape(u,M^2,1);

% Set initial guess
x0=zeros(M^2,1);

% Set tolerance 
tol=1e-3;

% Set maximum nummber of Jacobi iterations per time step
maxit=10000;

%--------------------my added code-----------------------------------------
dropTol = 1e-3;
opts=struct('type','ict','droptol',dropTol,'michol','off');
L = ichol(A,opts);
%--------------------------------------------------------------------------
% Loop over the timesteps
for n=1:count
    % Pause to give graphics routine time to catch up
    % pause(0.01);    
    % Solve linear system using pcg to advance one time step
    [x, flag, relres, it, resvec]=pcg(A,b,tol,maxit,L,L',x0);
    
    % Update right hand side
    b=x;
    
    % Display solution
    clf; u=reshape(x,M,M); u1=u(1:nb:end,1:nb:end);  surf(tau,tau,u1); ...
         view([0,90]); axis([0 1 0 1 0 maxtemp]);     
end
toc



