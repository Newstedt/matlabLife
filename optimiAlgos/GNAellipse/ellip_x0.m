function [c,ax,delta,theta]=ellip_x0(x,y,c,delta)
%ELLIP_X0 Calculate starting approximation for the ellipse problem.
%
%[c,ax,delta,theta]=ellip_x0(x,y[,c][,delta])
%
%Constructs and solves the linear problem
%x'*U*x+b'*x-gamma=0
%for U,b,gamma and converts to the geometric parameters cx,cy,...

% v1.0  2000-05-14. Niclas Borlin, niclas@cs.umu.se.
% v1.1  2000-05-19. Added options to specify c and/or delta.

switch (nargin)
case 2
	% x,y only

	% Right-hand side = gamma.
	g=ones(size(x));

	% Center around the mean.
	xm=mean(x);
	ym=mean(y);
	x=x-xm;
	y=y-ym;

	% Create and solve linear system.
	A=[x(:).^2,2*x(:).*y(:),y(:).^2,x(:),y(:)];
	v=A\g;

	% Unpack solution.
	U=[v(1:2),v(2:3)];
	b=v(4:5);

	% Convert to geometric parameters.
	[c,ax,delta]=mat2geo(U,b,1);

	% Use standard Q matrix.
	[dummy,dummy,dummy,Q]=geo2mat(c,ax,delta);
	
	% Calculate starting approximations of theta.
	csTheta=diag(1./ax)*Q'*[x(:)'-c(1);y(:)'-c(2)];
	theta=atan2(csTheta(2,:),csTheta(1,:));

	% Decenter from the origin.
	c=c+[xm;ym];
case 3
	if (length(c)==1)
		% x,y,delta
		delta=c;

		% Create rotation matrix.
		[dummy,dummy,dummy,Q]=geo2mat(zeros(2,1),ones(2,1),delta);

		% Center around the mean.
		xm=mean(x);
		ym=mean(y);
		x=x-xm;
		y=y-ym;
		xCentered=x;
		yCentered=y;

		% Rotate all points.
		xyRotated=Q'*[x(:)';y(:)'];

		x=xyRotated(1,:)';
		y=xyRotated(2,:)';

		% Create and solve linear system.
		A=[x.^2,y.^2,x,y];
		g=ones(size(x));
		v=A\g;
		
		% Unpack solution.
		L=diag(v(1:2));
		d=v(3:4);

		U=Q*L*Q';
		b=Q*d;
		
		[c,ax,dummy]=mat2geo(U,b,1);
		
		% Calculate starting approximations of theta.
		csTheta=diag(1./ax)*Q'*[xCentered'-c(1);yCentered'-c(2)];
		theta=atan2(csTheta(2,:),csTheta(1,:));

		c=c+[xm;ym];
	else
		% x,y,c

		% Create and solve linear system.
		A=[x.^2-2*c(1)*x,...
		   2*x.*y-2*c(1)*y-2*c(2)*x,...
		   y.^2-2*c(2)*y];
		g=ones(size(x));
		v=A\g;
		
		% Unpack solution.
		U=[v(1:2),v(2:3)];
		b=-2*U*c;
	
		% Convert to geometric parameters.
		[dummy,ax,delta]=mat2geo(U,b,1);
		
		[dummy,dummy,dummy,Q]=geo2mat(c,ax,delta);

		% Calculate starting approximations of theta.
		csTheta=diag(1./ax)*Q'*[x(:)'-c(1);y(:)'-c(2)];
		theta=atan2(csTheta(2,:),csTheta(1,:));
	end
case 4
	% x,y,c,delta

	% Create rotation matrix.
	[dummy,dummy,dummy,Q]=geo2mat(zeros(2,1),ones(2,1),delta);

	xOrig=x;
	yOrig=y;
	
	% Rotate all points.
	xyRotated=Q'*[x(:)';y(:)'];

	x=xyRotated(1,:)';
	y=xyRotated(2,:)';

	% Rotate center.
	d=Q'*c;
	
	% Create and solve linear system.
	A=[x(:).^2-2*d(1)*x,y(:).^2-2*d(2)*y];
	g=ones(size(x));
	v=A\g;
		
	% Unpack solution.
	L=diag(v(1:2));
	
	U=Q*L*Q';
	b=-2*U*c;
	[dummy_c,ax,dummy_delta]=mat2geo(U,b,1);	

	% Calculate starting approximations of theta.
	csTheta=diag(1./ax)*Q'*[xOrig(:)'-c(1);yOrig(:)'-c(2)];
	theta=atan2(csTheta(2,:),csTheta(1,:));
otherwise
	error('Illegal number of parameters')
end
