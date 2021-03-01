function [h1,h2,h3]=plotellipse_with_res(ax,x,Q)
%PLOTELLIPSE_WITH_RES Plot ellipse with residuals.
%
%   PLOTELLIPSE_WITH_RES(AX,X,Q) plots an ellipse with parameters from
%   X=[CX,CY,A,B,DELTA,TH1,...,THN]' in the axes AX. Three types of
%   lines are drawn: The full ellipse in blue, the estimated points
%   with red crosses, and the residual vectors in dashed green.
%
%Example:
%   pts=ellipse_data(1,1);
%   [c,ax,delta,theta]=ellipse_x0(pts(1,:)',pts(2,:)');
%   x0=[c(:);ax(:);delta;theta(:)];
%   tol=1e-6;
%   maxIter=10000;
%   d0=norm(x0);
%   mu=0.25;
%   eta=0.75;
%   [x,code,niter,X,deltas,rhos,steps]=levenberg('ellipse_res',x0,tol,maxIter,d0,mu,eta,{pts});
%   clf
%   n=size(X,2);
%   iStep=floor(sqrt(n));
%   ii=unique([1:iStep:n,max(n-iStep,1):n]);
%   for i=ii
%       plotellipse_with_res(gca,X(:,i),pts);
%       axis equal
%       title(sprintf('Iteration %d',i-1))
%       pause % or pause(0.01)
%   end

if ~ishandle(ax), error('First argument must be handle to an axes.'); end

if any(size(Q)~=[2,length(x)-5])
    error('Size error.');
end

% Unpack parameter vector.
c=x(1:2);
a=x(3);
b=x(4);
delta=x(5);
theta=x(6:end);
n=length(theta);

R=[cos(delta),-sin(delta);
   sin(delta),cos(delta)];
A=diag([a,b]);

% Points on full ellipse.
t=linspace(0,2*pi,101);
xyEll=repmat(c,1,length(t))+R*A*[cos(t);sin(t)];

% Add the ellipse line in blue.
h1=line(xyEll(1,:),xyEll(2,:),'parent',ax,'color','b','marker','none',...
        'linestyle','-');

% Estimated observations.
ct=cos(theta(:).');
st=sin(theta(:).');

xyEst=repmat(c,1,n)+R*A*[ct;st];

% Are estimated points already plotted?
estTag='estPts';
h2=findobj(ax,'type','line','tag',estTag);
if isempty(h2)
    % No, create one.
    h2=line(xyEll(1,:),xyEll(2,:),'parent',ax,'color','r','marker','x',...
            'linestyle','none','tag',estTag);
else
    % Yes, replace x/y-data.
    set(h2,'xdata',xyEll(1,:),'ydata',xyEll(2,:));
end

% Residuals
xRes=reshape([xyEst(1,:);Q(1,:);nan(1,n)],[],1);
yRes=reshape([xyEst(2,:);Q(2,:);nan(1,n)],[],1);

% Does residual lines already exist?
resTag='residuals';
h3=findobj(ax,'type','line','tag',resTag);
if isempty(h3)
    % No, create one.
    h3=line(xRes,yRes,'parent',ax,'color',0.5*[0,1,0],'marker','none',...
            'linestyle','--','tag',resTag);
else
    % Yes, replace x/y-data.
    set(h3,'xdata',xRes,'ydata',yRes);
end