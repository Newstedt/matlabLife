function [xy,dIO,dEO,dOP]=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%PM_MULTIEULERPINHOLE1 Pinhole camera projection into multiple camera stations.
%
%[xy,dIO,dEO,dOP]=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis[,cIO,cEO,cOP])
%IO     - matrix with camera inner orientation as columns
%         pp - principal point [xp;yp].
%         f  - focal length.
%         K  - radial lens distortion coefficients.
%         P  - tangential lens distortion coefficients.
%         a  - with affine lens distortion coefficients.
%         u  - image unit.
%nK     - number of radial koefficients.
%nP     - number of tangential koefficients.
%EO     - matrix with photo camera outer orientation as columns
%         C     - camera center [Xc,Yc,Zc].
%         ang   - Euler angles.
%         ax    - Euler axis sequence
%                 0 - 'xyz'.
%                 1 - 'zxz'.
%cams   - vector of camera numbers for each photo.
%OP     - 3 x nObj matrix with object coordinates
%vis    - sparse matrix indicating if obj point i is visible in photo j.
%cIO    - should we calculate partial derivatives w.r.t. internal
%         orientation? Scalar or matrix of same size as IO.
%cEO    - should we calculate partial derivatives w.r.t. external
%         orientation? Scalar or matrix of same size as EO.
%cOP    - should we calculate partial derivatives w.r.t. object points?
%xy     - projected points [x;y](:).
%dIO    - jacobian w.r.t. internal orientation.
%dEO    - jacobian w.r.t. external orientation.
%dOP    - jacobian w.r.t. object points.

% v1.0  2004-04-14. Niclas Borlin, niclas@cs.umu.se.
%                   Adapted from pm_multieulercam1.m.

if (nargin<8), cIO=(nargout>1); end
if (nargin<9), cEO=(nargout>2); end
if (nargin<10), cOP=(nargout>3); end

dIO=[];
dEO=[];
dOP=[];

nCam=size(IO,2);
nPhotos=size(EO,2);
nObjs=size(OP,2);

if (length(cams)==1)
	% Same camera for all photos.
	cams=repmat(cams,nPhotos,1);
end

% Total number of projected points.
nProj=nnz(vis);

if (all(~cIO(:)) & all(~cEO(:)) & all(~cOP(:)))
	% No partial derivatives at all.
	
	% Preallocate point matrix for speed.
	xy=zeros(2*nProj,1);
	xyBase=1;

	for i=1:nPhotos
		% Get camera station.
		camStation=EO(:,i);
		center=camStation(1:3);
		ang=camStation(4:6);
		if (camStation(7)==0)
			seq='xyz';
		else
			seq='zxz';
		end

		% Get inner orientation.
		camNo=cams(i);
		[pp,f,K,P,a,u]=UnpackIO(IO(:,camNo),nK,nP);
	
		% Get object points visible in this image
		v=vis(:,i);
		obj=OP(:,v);
	
		% Project into image.
		imPt=pm_eulerpinhole1(pp,f,obj,center,ang,seq);
	
		% Find out where to store points.
		[ixPt,xyBase]=pindex(nnz(v)*2,xyBase);
		xy(ixPt)=imPt(:);
	end
else
	% Which IO partial derivatives are requested?
	if (length(cIO)==1)
		cIO=repmat(cIO,size(IO,1),size(IO,2));
	end

	% Preallocate IO jacobian.
	% Number of wanted parameters
	ioCols=nnz(cIO);
	% Max number of non-zero elements.
	ioMaxNnz=nProj*2*max(sum(cIO));
	dIO=sparse([],[],[],nProj*2,ioCols,ioMaxNnz);

	% Create arrays of columns indices for IO derivatives.
	[ixpp,ixf,ixK,ixP,ixa,ixu]=CreateIOColumnIndices(cIO,nK,nP);
	
	
	% Which EO partial derivatives are requested?
	if (length(cEO)==1)
		cEO=repmat(cEO,6,nPhotos);
	end

	% Preallocate EO jacobian.
	% Number of wanted parameters
	eoCols=nnz(cEO);
	% Max number of non-zero elements.
	eoMaxNnz=nProj*2*max(sum(cEO));
	dEO=sparse([],[],[],nProj*2,eoCols,eoMaxNnz);

	% Create arrays of columns indices for EO derivatives.
	[ixC,ixAng]=CreateEOColumnIndices(cEO);

	
	% Which OP partial derivatives are requested?
	if (length(cOP)==1)
		cOP=repmat(cOP,3,nObjs);
	end

	% Preallocate OP jacobian.
	% Number of wanted parameters
	opCols=nnz(cOP);
	% Max number of non-zero elements.
	opMaxNnz=nProj*2*max(sum(cOP));
	if (0)
		dOP=sparse([],[],[],nProj*2,opCols,opMaxNnz);
	else
		dOP=zeros(nProj*2,opCols);
	end

	% Create array of columns indices for OP derivatives.
	ixOP=CreateOPColumnIndices(cOP);
	
	% Preallocate point matrix for speed.
	xy=zeros(2*nProj,1);
	xyBase=1;

	for i=1:nPhotos
		% Get camera station.
		camStation=EO(:,i);
		center=camStation(1:3);
		ang=camStation(4:6);
		if (camStation(7)==0)
			seq='xyz';
		else
			seq='zxz';
		end

		% Which outer orientation parameters are interesting?
		cC=cEO(1:3,i);
		cAng=cEO(4:6,i);
		
		% Get inner orientation.
		camNo=cams(i);
		[pp,f,K,P,a,u]=UnpackIO(IO(:,camNo),nK,nP);
	
		% Which inner orientation parameters are interesting?
		[cpp,cf,cK,cP,ca,cu]=UnpackIO(cIO(:,camNo),nK,nP);
		
		% Get object points visible in this image
		v=vis(:,i);
		obj=OP(:,v);
		% ...and which partial derivatives to calculate...
		cObj=cOP(:,v);
		% ... and in which columns to store them.
		ixObj=ixOP(:,v);
		
		% Project into image.
		[imPt,dpp,df,dO,dC,dAng]=...
			pm_eulerpinhole1(pp,f,obj,center,ang,seq,...
						 any(cpp),cf,cObj,any(cC),any(cAng));
	
		% Find out where to store points.
		[ixPt,xyBase]=pindex(nnz(v)*2,xyBase);
		xy(ixPt)=imPt(:);
		
		% IO jacobians.
		if (any(cpp))
			dIO(ixPt,ixpp(cpp,camNo))=dpp(:,cpp);
		end
		if (cf)
			dIO(ixPt,ixf(camNo))=df;
		end
		
		% EO jacobians.
		if (any(cC))
			dEO(ixPt,ixC(cC,i))=dC(:,cC);
		end
		if (any(cAng))
			dEO(ixPt,ixAng(cAng,i))=dAng(:,cAng);
		end

		% OP jacobians.
		if (any(cObj(:)))
			dOP(ixPt,ixObj(cObj))=dO;
		end
	end

	% Verify
	%disp([nnz(dIO),ioMaxNnz]);
	%disp([nnz(dEO),eoMaxNnz]);
end


function [pp,f,K,P,a,u]=UnpackIO(IO,nK,nP)
%Unpack inner orientation vector.

pp=IO(1:2);
f=IO(3);
K=IO(3+[1:nK]);
P=IO(3+nK+[1:nP]);
a=IO(3+nK+nP+[1:2]);
u=IO(3+nK+nP+2+1);

function [ixpp,ixf,ixK,ixP,ixa,ixu]=CreateIOColumnIndices(cIO,nK,nP);
% Create arrays of columns indices for IO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

% How many cameras do we have?
nCams=size(cIO,2);

ix=reshape(cumsum(cIO(:)),6+nK+nP,nCams).*cIO;

ixpp=ix(1:2,:);
ixf=ix(3,:);
ixK=ix(3+[1:nK],:);
ixP=ix(3+nK+[1:nP],:);
ixa=ix(3+nK+nP+[1:2],:);
ixu=ix(3+nK+nP+2+1,:);

function [ixC,ixAng]=CreateEOColumnIndices(cEO);
% Create arrays of columns indices for EO derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

if (size(cEO,1)>6)
	cEO=cEO(1:6,:);
end

% How many photos do we have?
nPhotos=size(cEO,2);

ix=reshape(cumsum(cEO(:)),6,nPhotos).*cEO;

ixC=ix(1:3,:);
ixAng=ix(4:6,:);

function ixOP=CreateOPColumnIndices(cOP);
% Create arrays of columns indices for OP derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

% How many points do we have?
nObjs=size(cOP,2);

ixOP=reshape(cumsum(cOP(:)),3,nObjs).*cOP;
