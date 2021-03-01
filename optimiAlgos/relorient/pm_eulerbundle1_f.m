function [f,J,JJ]=pm_eulerbundle1_f(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%PM_EULERBUNDLE1_F Residual function for Euler camera bundle adjustment.
%
%[f,J,JJ]=pm_eulerbundle1_f(x,pt,ptCams,IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP)
%x      - vector of unknowns.
%pt     - measured projected points.
%ptCams - physical camera for each point.
%IO     - vector with camera inner orientation as columns
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
%cIO    - logical matrix indicating which elements of IO are unknown.
%cEO    - logical matrix indicating which elements of EO are unknown.
%cOP    - logical matrix indicating which elements of OP are unknown.
%f      - residual vector.
%J      - jacobian w.r.t. the unknowns in order [dEO,dOP]
%JJ     - numerical approximation of J.

% v1.0  2004-04-12. Niclas Borlin, niclas@cs.umu.se.

% Calculate indices for unknown parameters.
base=1;

[ixIO,base]=pindex(nnz(cIO),base);
[ixEO,base]=pindex(nnz(cEO),base);

if (prod(size(cOP))==1)
	cOP=repmat(cOP,size(OP));
elseif (size(cOP,1)<size(OP,1))
	cOP=repmat(cOP,size(OP,1),1);
end
[ixOP,base]=pindex(nnz(cOP),base);

if (size(cEO,1)<size(EO,1))
	cEO(end+1,1)=0;
end

% Copy the current approximations of the unknown values.
IO(cIO)=x(ixIO);
EO(cEO)=x(ixEO);
OP(cOP)=x(ixOP);

if (nargout>2)
	JJ=jacapprox(mfilename,x,1e-6,{pt,ptCams,IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP});
end

if (nargout<2)
	% Project into pinhole camera.
	xy=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis);

	% Remove lens distortion from measured points.
	ptCorr=pm_multilenscorr1(pt,IO,nK,nP,ptCams,size(IO,2));
	
	f=xy(:)-ptCorr(:);
else
	% Project into pinhole camera.
	[xy,dIO1,dEO,dOP]=pm_multieulerpinhole1(IO,nK,nP,EO,cams,OP,vis,cIO,cEO,cOP);
	
	% Remove lens distortion from measured points.
	[ptCorr,dIO2]=pm_multilenscorr1(pt,IO,nK,nP,ptCams,size(IO,2),cIO);

	f=xy(:)-ptCorr(:);
	
	J=[dIO1-dIO2,dEO,dOP];
end

