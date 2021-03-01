function [xy,dIO,dp]=pm_multilenscorr1(p,IO,nK,nP,cams,nCams,cIO,cp)
%PM_MULTILENSCORR1 Correct for lens distortion in multiple images.
%
%[xy,dIO,dp]=pm_multilenscorr1(p,IO,nK,nP,cams,nCams[,cIO,cp])
%p      - 2 x n matrix with measured points.
%IO     - matrix with camera inner orientation as columns
%         pp - principal point [xp;yp].
%         f  - focal length.
%         K  - radial lens distortion coefficients.
%         P  - tangential lens distortion coefficients.
%         a  - with affine lens distortion coefficients.
%         u  - image unit.
%nK     - number of radial koefficients.
%nP     - number of tangential koefficients.
%cams   - vector of camera numbers for each point.
%nCams  - total number of cameras.
%cIO    - should we calculate partial derivatives w.r.t. internal
%         orientation? Scalar or matrix of same size as IO.
%cp     - should we calculate partial derivatives w.r.t. points?
%         2 x n or 1 x n vector.
%xy     - points corrected for lens distortion [x;y](:).
%dIO    - jacobian w.r.t. internal orientation.
%dp     - jacobian w.r.t. 2d points.

% v1.0  2004-04-10. Niclas Borlin, niclas@cs.umu.se.

if (nargin<7), cIO=(nargout>1); end
if (nargin<8), cp=(nargout>2); end

% Total number of points.
nPts=size(p,2);

dp=[];
dIO=sparse([],[],[],nPts*2,nnz(cIO));

if (length(cams)==1)
	% Same camera for all points.
	cams=repmat(cams,nPts,1);
end

if (all(~cIO(:)) & ~(any(cp)))
	% No partial derivatives.
	
	% Preallocate point matrix for speed.
	xy=zeros(size(p));

	for i=1:nCams
		% Get inner orientation.
		[pp,f,K,P,a,u]=UnpackIO(IO(:,i),nK,nP);
	
		% Get points taken with this camera.
		ix=cams==i;
		
		% Convert to mm.
		q=p(:,ix)/u;
		
		% Lens distortion.
		lens=pm_lens1(q,pp,K,P);
		
		% Correct for lens distortion.
		xy(:,ix)=q-lens;
	end
	xy=xy(:);
else
	% Which partial derivatives are requested?
	if (length(cIO)==1)
		cIO=repmat(cIO,size(IO,1),size(IO,2));
	end

	if (prod(size(cp))==1)
		cp=repmat(cp,1,nPts);
	end
	
	% Preallocate jacobians.
	
	% Number of wanted internal parameters
	ioCols=nnz(cIO);
	% Max number of non-zero elements.
	ioMaxNnz=nPts*2*max(sum(cIO));
	dIO=sparse([],[],[],nPts*2,ioCols,ioMaxNnz);

	% Create arrays of columns indices for IO derivatives.
	[ixpp,ixf,ixK,ixP,ixa,ixu]=CreateIOColumnIndices(cIO,nK,nP);

	% Number of wanted points.
	ptCols=2*nnz(cp);
	ptMaxNnz=4*nPts;
	dp=sparse([],[],[],nPts*2,ptCols,ptMaxNnz);
	
	% Create array of columns indices for point derivatives.
	ixCol=CreatePtColumnIndices(cp);

	% Preallocate point matrix for speed.
	xy=zeros(size(p));

	for i=1:nCams
		% Get inner orientation.
		[pp,f,K,P,a,u]=UnpackIO(IO(:,i),nK,nP);
	
		% Which inner orientation parameters are interesting?
		[cpp,cf,cK,cP,ca,cu]=UnpackIO(cIO(:,i),nK,nP);
		
		% Get points taken with this camera.
		ix=find(cams==i);

		% Convert to mm.
		q=p(:,ix)/u;
		
		if (any(cp(:)))
			dqdp=speye(prod(size(q)))./u;
		else
			dqdp=[];
		end
		if (cu)
			dqdu=reshape(-p(:,ix)./u^2,prod(size(q)),1);
		else
			dqdu=[];
		end
		
		% Lens distortion.
		calcID=cu | any(cp(:));
		[lens,dldq,dldpp,dldK,dldP]=pm_lens1(q,pp,K,P,calcID,...
											 any(cpp),any(cK),any(cP));
		
		% Correct for lens distortion.
		xy(:,ix)=q-lens;

		% Calculate row indices in jacobian.
		ixRow=[(ix(:)-1)*2+1,(ix(:)-1)*2+2]';
		ixRow=ixRow(:);
		
		% IO jacobians.
		if (any(cpp))
			dIO(ixRow,ixpp(cpp,i))=-dldpp(:,cpp);
		end
		if (cf)
			% Focal length does not take part in equation.
			% dIO(ixRow,ixf(i))=0;
		end
		if (any(cK))
			dIO(ixRow,ixK(cK,i))=-dldK(:,cK);
		end
		if (any(cP))
			dIO(ixRow,ixP(cP,i))=-dldP(:,cP);
		end
		if (any(ca))
			% Not implemented yet.
			% dIO(ixRow,ixa(ca,i))=0
		end
		if (cu)
			dIO(ixRow,ixu(i))=dqdu-dldq*dqdu;
		end
		if (any(cp(ix)))
			colix=ixCol(:,ix);
			dqdp=dqdp(:,colix~=0);
			dp(ixRow,colix(colix~=0))=dqdp-dldq*dqdp;
		end
	end

	% Verify
	%disp([nnz(dIO),ioMaxNnz]);
	%disp([nnz(dp),ptMaxNnz]);
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

function ixp=CreatePtColumnIndices(cp);
% Create arrays of columns indices for OP derivatives.
% A zero element means that the partial derivative should not be
% calculated/stored.

cp=[cp;cp];

% How many points do we have?
nPts=size(cp,2);

ixp=reshape(cumsum(cp(:)),2,nPts).*cp;
