format compact

% Change this for different camera pairs.
cpSelect=1;

% Put the name of your algorithm as third name.
algNames={'LS','C','LM'};
% Change last element to false if algorithm is trust-region.
algIsLineSearch=[true,true];

% Avoid reloading data.
if (exist('prob')==0)
    prob=loadpm('zurich.txt');
    pause(0.1);
end

% Test constants.
doSim=1;
doPlot=1;
doWrite=0;
keepFraction=1;
doWait=true;
use2ndaxes=1;
plotx0pts=0;
plotCameraTrace=1;
usePointHeuristics=1;
useDlt=0;

pairs=true;

disp(['Relative orientation example with Zurich City Hall data set'])

views=nan*ones(90,2);
views(1,:)=[45,8];
views(2,:)=[7,18];
views(3,:)=[4,26];
views(4,:)=[-20,24];
views(5,:)=[-11,26];
views(6,:)=[-8,24];
views(7,:)=[83,75];
views(8,:)=[10,10];
views(9,:)=[3,-8];
views(10,:)=[-1,14];
views(11,:)=[1,22];
views(12,:)=[0,22];
views(13,:)=[5,14];
views(14,:)=[4,12];
views(15,:)=[0,16];
views(16,:)=[-2,6];
views(17,:)=[0,18];
views(18,:)=[17,44];
views(19,:)=[-5,6];
views(20,:)=[-7,10];
views(21,:)=[-12,12];
views(23,:)=[75,90];
views(24,:)=[6,22];
views(25,:)=[5,22];
views(26,:)=[2,32];
views(27,:)=[26,90];
views(28,:)=[23,20];
views(30,:)=[70,-6];
views(31,:)=[-2,14];
views(32,:)=[48,30];
views(34,:)=[67,6];

axisSizes=nan*ones(90,6);
%axisSizes(1,:)=[-0.0266    0.9928   -0.1257    0.9644   -0.4890    0.1950];
axisSizes(6,:)=[-2.2862    1.7314   -0.9194    2.2493   -0.7310    0.8534];
axisSizes(18,:)=[-1.6341    1.2094   -0.7284    1.4043   -1.4757    1.3678];
axisSizes(22,:)=[-61.6737   46.4067  -27.1973   53.8630  -52.3000   55.7804];
axisSizes(25,:)=[-10.2101   10.3733   -0.4577   34.0494   -3.9777    5.9318];
axisSizes(32,:)=[-1.1952   1.5047   -0.3911    3.4409   -0.3475    0.4456];

rotData=zeros(90,3);
rotData(18,:)=[0,0,90];
rotData(22,:)=[0,0,90];

if (any([doSim,doPlot,doWrite,keepFraction,doWait,cpSelect,use2ndaxes,plotx0pts,plotCameraTrace,usePointHeuristics]))
end

x0Col=0.75*ones(1,3);

if (doSim)
    imSz=[1280,1024];
    if (exist(prob.cameras(1).imName))
	info=imfinfo(prob.cameras(1).imName);
	imSz=[info.Width,info.Height];
    end
    disp(sprintf('Assuming image size of %dx%d pixels',imSz));

    % Construct internal parameters.
    f=prob.cameras(1).inner(1);
    wh=prob.cameras(1).inner(4:5);
    pp=prob.cameras(1).inner(2:3);
    pp(2)=-pp(2);
    K=-prob.cameras(1).inner(6:7);
    P=-prob.cameras(1).inner(9:10);
    nK=length(K);
    nP=length(P);
    a=zeros(2,1);
    if (diff(imSz./wh)>1)
	warning('Different image ratio');
    end

    u=mean(imSz./wh);

    IO=[pp(:);f;K(:);P(:);a(:);u];

    camId=[4 5 12 13 14 15 24    25    26    27    28    29    30    31]';

    camPairs=[
        2     5   674
        2     4   667
        4     5   659
        1     8   294
        5     6   283
        11    14   263
        9    10   260
        6    11   259
        6    14   259
        1     9   257
        1     2   249
        5    11   245
        5    14   245
        1    13   211
        1    12   209
        2    12   205
        2    13   200
        7     8   176
        1     7   172
        8     9   167
        1    10   166
        12    13   165
        10    11   157
        2     3   121
        3     4   120
        3     5   111
        9    11    85
        4     6    65
        7     9    52
        2     6    51
        1    11    49
        3     6    43
        8    10    42
        1     5    33
        4    11    33
        6    10    32
        2    11    32
        2    14    32
        4    14    32
        10    14    32
        5    10    31
        5    12    31
        8    13    31
        2     8    28
        7    13    28
        1     4    27
        2     7    26
        4    12    26
        7    11     7
        1     6     4
        1     3     3
        3    11     3
        8    11     3
        1    14     3
        3    12     2
        3    14     2
        6    12     1
             ];

    camTriplets=[
        2     5     6   918
        1     2     5   892
        1     2     4   889
        2     5    11   889
        2     5    14   889
        4     5     6   883
        4     5    14   874
        4     5    11   873
        2     5    12   848
        2     4    12   846
        2     3     5   740
        2     4     5   732
        2     3     4   728
        3     4     5   726
        2     4     6   691
        2     4    11   668
        1     4     5   667
        4     5    12   666
        1     2     8   515
        1     8    13   474
        1     8    10   438
        1     8     9   430
        1     9    10   415
        1     7     9   399
        1     2     7   395
        10    11    14   388
        8     9    10   385
        6    10    11   384
        5    10    11   371
        1     2     3   367
        3     5     6   367
        1     7    13   355
        3     5    11   355
        9    10    11   342
        1     8    11   340
        1     9    11   329
        2     3    12   324
        1     7     8   318
        1     5     6   314
        1    11    14   309
        1     6    11   308
        7     8     9   307
        1     2     6   300
        3     6    14   300
        3     6    11   299
        5     6    11   297
        5     6    14   297
        4     6    14   292
        4     6    11   291
        1    10    11   290
        5     6    10   284
        2     6    11   278
        2     6    14   278
        1     2    13   264
        3    11    14   264
        4    11    14   264
        2    11    14   263
        1     6    14   262
        1    12    13   259
        1     2    12   257
        2    12    13   252
        5    10    14   246
        1     7    11   220
        2     8    13   203
        2     7    13   202
        1     6    10   198
        1    10    14   197
        7     8    11   182
        7     8    13   179
        2     7     8   178
        2     3     6   173
        3     4     6   158
        2     3    11   152
        2     3    14   151
        3     4    11   150
        3     4    14   150
        1     3     5   143
        7     9    11   134];

    % Keep camera pairs with >=42 common points.
    keep=camPairs(:,3)>=42;
    %camPairs=camPairs(keep,:);

    %cpSelect=floor(rand*size(camPairs,1))+1;
    %cpSelect=cpSelect+1
    if (pairs)
        useCams=camPairs(cpSelect,1:2);
        disp(sprintf('Cameras %d %d (camIds %d %d)',useCams,camId(useCams)));
    else
        useCams=camTriplets(cpSelect,1:3);
        disp(sprintf('Cameras %d %d %d (camIds %d %d %d)',useCams,camId(useCams)));
    end    

    %useCams=uc(1:2)

    %useCams=9:10

    camSeq=[14,5,13,15,28,31,4,25,26,27,30,29,24,12];
    %useCams=find(ismember(camId,camSeq(1:4)));

    %useCams=1:length(prob.cameras);

    % Data storage:
    %
    % IO internal parameters.
    % EO external parameters, columnwise.
    %
    % OPid contains object point ids
    % OP   contains object point XYZ, sorted by id.
    %
    % isCtrl shows if OP is control point
    %
    % vis(i,j)=1 if object point i is visible in image j
    %
    % markPos mark positions sorted by image, id.

    % Camera stations.
    cams=cat(1,prob.cameras.outer);
    C=cams(:,1:3)';
    ang=cams(:,[6,5,4])'/180*pi;

    realC=C;

    % Construct exterior orientation parameters.
    EO=[C;ang];
    % Indicate all rotations are Euler 'xyz'.
    EO(7,1)=0;

    EOeuler=EO;

    % Camera calibration matrix.
    K=diag([-f,-f,1]);
    K(1:2,3)=pp';

    realCamMatrices=repmat(nan,12,size(EO,2));
    for i=1:size(EO,2)
        Ci=EO(1:3,i);
        Ri=rotmat(EO(4:6,i));
        realCamMatrices(:,i)=reshape(K*Ri*[eye(3),-Ci],12,1);
    end
    
    % Construct exterior orientation parameters.
    EOrq=[C;zeros(9,size(ang,2))];
    for i=1:size(EO,2)
	EOrq(4:12,i)=unroll(rotmat(ang(:,i),'xyz'));
    end

    % Object points.
    OP=prob.objPts(:,2:4)';
    OPid=prob.objPts(:,1);

    % Sort according to id.
    [dummy,ix]=sort(OPid);
    OPid=OPid(ix);
    OP=OP(:,ix);
    if (isempty(prob.ctrlPts))
	isCtrl=logical(zeros(size(OPid)));
    else
	isCtrl=ismember(OPid,prob.ctrlPts(:,1));
    end

    realOP=OP;

    % Sort mark points after image, then id.
    markPts=msort(prob.markPts);

    % Remove mark points in unused images.
    imNo=markPts(:,1)'+1;
    unused=~ismember(imNo,useCams);
    markPts(unused,:)=[];

    % Remove unreconstructed points.
    markId=markPts(:,2);
    reconstructed=ismember(markId,OPid);
    markPts=markPts(reconstructed,:);

    % Remove points visible in too few images.
    imNo=markPts(:,1)'+1;
    markId=markPts(:,2);
    vis=logical(sparse(markId,imNo,1)); % Mapping from (id,image)
    if (size(vis,1)<max(OPid))
	vis(max(OPid),1)=0;
    end
    if (size(vis,2)<max(useCams))
	vis(1,max(useCams))=0;
    end

    vis=vis(OPid,useCams); % Mapping from (OPid(i),image(j))
    ptOK=sum(vis')>1;

    nKeep=8+floor(keepFraction*(nnz(ptOK)-8));

    % Pick a selection of the OK points.
    zz=find(ptOK);
    zzi=zz(randperm(length(zz)));

    ptOK(zzi(nKeep+1:end))=0;

    OPid=OPid(ptOK);
    OP=OP(:,ptOK);
    isCtrl=isCtrl(ptOK);

    markId=markPts(:,2);
    markPtOK=ismember(markId,OPid);
    markPts=markPts(markPtOK,:);

    imNo=markPts(:,1)'+1;
    markId=markPts(:,2);
    markPos=markPts(:,3:4)';

    % Construct visibility matrix.
    vis=sparse(markId,imNo,1); % Mapping from (id,image)
    if (size(vis,1)<max(OPid))
	vis(max(OPid),1)=0;
    end
    if (size(vis,2)<max(useCams))
	vis(1,max(useCams))=0;
    end
    vis=logical(vis(OPid,useCams)); % Mapping from (OPid(i),image(j))

    realOP=OP;
    realVis=vis;

    ptAngles=zeros(size(realOP,2),1);
    ptDist=zeros(size(realOP,2),2);
    for i=1:size(realOP,2)
	v=realOP(:,i);
	c1=realC(:,useCams(1));
	c2=realC(:,useCams(2));
	ptAngles(i)=rad2deg(subspace([c1-v],[c2-v]));
	ptDist(i,1)=norm(v-c1);
	ptDist(i,2)=norm(v-c2);
    end
    
    sv=full(unique(sum(vis')));
    for i=1:length(sv)
        disp(sprintf('Number of points seen in %d images: %d',...
                     sv(i),nnz(sum(vis')==sv(i))));
    end

    % Create mapping from OPid,image to column in markPts.
    colPos=reshape(cumsum(vis(:)),size(vis)).*vis;

    % Number of object points.
    n=size(OP,2);

    % Number of cameras.
    m=size(EO,2);

    disp('Using nominal -y')
    markPos(2,:)=-markPos(2,:);

    % Remove lens distortion and convert to mm.
    correctedPt=roll(pm_multilenscorr1(markPos,IO,nK,nP,1,1),2);

    disp('Using real points');

    disp('Initializing')
    OPsave=OP;
    EOsave=EO;
    IOsave=IO;

    % Determine camera orientations.
    XForm=rotmat(deg2rad([-90,0,0]));

    allCommon=all(vis');

    % Set up object point matrix.
    OP=repmat(nan,size(OP));
    % Are the reconstructed object points ok?
    OPok=repmat(false,1,size(OP,2));
    % How many cameras were used to reconstruct each object point?
    OPcams=zeros(1,size(OP,2));

    % Construct a camera-point coincidence matrix.
    nCommon=sparse(length(useCams),length(useCams));
    for i=1:length(useCams)
        for j=i+1:length(useCams)
            % Number of points visible in camera i and j.
            nCommon(i,j)=nnz(vis(:,i) & vis(:,j));
            nCommon(j,i)=nnz(vis(:,i) & vis(:,j));
        end
    end

    % Start with one of the cameras with the most common points.
    oriented=repmat(false,1,length(useCams));
    [dummy,r,c]=max2(nCommon);
    if (sum(nCommon(r,:))>sum(nCommon(c,:)))
        i=r;
    else
        i=c;
    end

    camMatrices=repmat(nan,12,length(useCams));
    % Initialize reference camera.
    oriented(i)=true;
    camMatrices(:,i)=reshape(eye(3,4),12,1);

    EOnew=repmat(nan,7,length(useCams));
    EOnew(end,:)=0;
    EOnew(1:6,i)=0;
    baseCam=i;

    inFrontOf=vis;

    % Repeat until all camera have been oriented.
    while (~all(oriented))
        % Decide if we should use relative orientation or DLT.
        
        % Find the unoriented camera with the largest number of points common
        % with an oriented camera.
        ixOriented=find(oriented);
        ixUnoriented=find(~oriented);
        
        [relOrientPts,r,c]=max2(nCommon(ixOriented,ixUnoriented));

        % Find the unoriented camera with largest number of visible object
        % points.
        [dltPoints,i]=max(sum(vis(OPok,~oriented)));
        
        if (~useDlt | dltPoints==0)
            refCam=ixOriented(r);
            newCam=ixUnoriented(c);
            %       camId(useCams([refCam,newCam]))'
            
            % Points visible in both images.
            ptsInPair=vis(:,refCam) & vis(:,newCam);
            % Points visible in both images and with previous 3d data.
            common=ptsInPair' & OPok;

            % Determine fundamental matrix from points.
            pt1=correctedPt(:,colPos(ptsInPair,refCam));
            pt2=correctedPt(:,colPos(ptsInPair,newCam));
            if (1)
                [P1,P2,X,inFront]=pm_relorient(homogenous(pt1),homogenous(pt2),K,K,false);
            else
                pt1n=K\homogenous(pt1);
                pt2n=K\homogenous(pt2);
                Evec=essmat5(pt2n,pt1n);
                bestSp=-inf;
                bestErr=inf;
                for i=1:size(Evec,2)
                    [P1new,P2new,Xnew,inFrontnew,sp]=camsfrome(reshape(Evec(:,i),3,3),pt1n,pt2n);
                    if (sp(1)-sp(2)>bestSp)
                        bestSp=sp(1)-sp(2);
                        P1=P1new;
                        P2=P2new;
                        X=Xnew;
                        inFront=inFrontnew;
                        bestErr=sum(sum((euclidian(P1*X)-euclidian(pt1n)).^2))+...
                                sum(sum((euclidian(P2*X)-euclidian(pt2n)).^2));
                    elseif (sp(1)-sp(2)==bestSp)
                        newErr=sum(sum((euclidian(P1new*Xnew)-euclidian(pt1n)).^2))+...
                               sum(sum((euclidian(P2new*Xnew)-euclidian(pt2n)).^2));
                        if (newErr<bestErr)
                            P1=P1new;
                            P2=P2new;
                            X=Xnew;
                            inFront=inFrontnew;
                            bestErr=newErr;
                        end
                    end
                end
            end
            X=normhomo(X);
            
            % Determine scale of object points already known.
            if (nnz(oriented)<2)
                scale=1;
                Tscale=eye(4);
            elseif (nnz(common)<2)
                % Too few common points.
                error('Too few common points');
            else
                oldDist=disttbl(OP(:,common));
                scaleOld=sqrt(mean(oldDist(logical(triu(ones(size(oldDist)),1))).^2));
                newDist=disttbl(X(1:3,common(ptsInPair)));
                scaleNew=sqrt(mean(newDist(logical(triu(ones(size(newDist)),1))).^2));
                scale=scaleOld/scaleNew;
                Tscale=diag([1,1,1,1/scale]);
            end
            
            % Scale points and camera 2.
            Xs=Tscale*X;
            P2s=P2*inv(Tscale);
            
            % Create transformation that would transform reference camera P1 to
            % its real position.
            P1r=reshape(camMatrices(:,refCam),3,4);
            T=[P1(:,1:3)\P1r(:,1:3), P1(:,1:3)\(P1r(:,4)-P1(:,4));0,0,0,1];

            % Check that T does it's job.
            %norm(P1r-P1*T,'fro')
            
            % Apply transformation on P2 and inverse transformation to points.
            P2r=P2s*T;
            Xr=normhomo(T\Xs);
            
            OPok(ptsInPair)=OPok(ptsInPair) | inFront;
            OP(:,ptsInPair)=Xr(1:3,:);
        else
            newCam=ixUnoriented(i);
            % What object points can we trust?
            OPix=vis(:,newCam) & OPok';
            % Estimate camera matrix with the DLT.
            P2=dlt(K\homogenous(correctedPt(:,colPos(OPix,newCam))),OP(:,OPix));
            % Unpack camera.
            [lambda,P2c,M,pp,f]=melendecomp(P2);
            [rad2deg(derotmat3d(M)),det(M)]
            [U,S,V]=svd(P2(:,1:3));
            M2=U*V'*det(U*V');
            [rad2deg(derotmat3d(M2)),det(M2)]
            MM=K\reshape(realCamMatrices(:,useCams(i)),3,4);
            [rad2deg(derotmat3d(MM(:,1:3))),det(MM(:,1:3))]
            P2r=M2*[eye(3),-P2c];
        end
        
        % Store camera matrix.
        camMatrices(:,newCam)=reshape(P2r,12,1);
        
        EOnew(1:3,newCam)=euclidian(null(P2r));
        EOnew(4:6,newCam)=derotmat3d(P2r(:,1:3))';
        
        % Mark new camera as oriented.
        oriented(newCam)=true;

        for j=find(sum(vis(:,oriented),2)>1)'
            PP=[];
            xy=[];
            for k=find(vis(j,:) & oriented)
                PP=[PP;K*reshape(camMatrices(:,k),3,4)];
                xy=[xy,correctedPt(:,colPos(j,k))];
            end
            pt=pm_forwintersect1(PP,xy,'euclidian');
            OP(:,j)=pt;
            OPok(j)=true;
        end

        % Find out which points are in front of which cameras.
        inFrontOf=vis;
        for i=find(oriented)
            P=K*reshape(camMatrices(:,i),3,4);
            inFront=ptdepth(P,OP(:,vis(:,i)))<0;
            inFrontOf(vis(:,i),i)=inFront;
            OPok(vis(:,i))=OPok(vis(:,i)) & inFront;
        end
        
        if (0)
            set(figure(1),'position',[0,0,560,420])
            plot3(OP(1,OPok),OP(2,OPok),OP(3,OPok),'b.');
            line(OP(1,~OPok),OP(2,~OPok),OP(3,~OPok),'color','r','marker','.','linestyle','none');
            line(EOnew(1,:),EOnew(2,:),EOnew(3,:),'marker','o','linestyle', ...
                 'none','color','b');
            hold on
            for i=find(oriented)
                camSize=1;
                [cam,camCol]=cameraicon(camSize*[0.11,0.08,0.04]);
                cam(:,:,3)=-cam(:,:,3);
                [m,n,p]=size(cam);
                RR=pm_eulerrotmat(EOnew(4:6,i));
                CC=EOnew(1:3,i);
                T=[RR*[eye(3),-CC];0,0,0,1];
                cam=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);
                surf(cam(:,:,1),cam(:,:,2),cam(:,:,3),camCol,'edgealpha',1, ...
                     'facealpha',0.15);
                colormap([1,0,0;0,1,0;0.5,0.5,1]);
            end
            hold off
            axis equal
            axis vis3d
            %pause
            drawnow
        end
    end

    camMatricesSave=camMatrices;

    EOeuler=EOnew;

    cOP=logical(ones(size(OP)));
    %cOP(:,isCtrl)=0;

    referenceBase=max(abs(EOeuler(1:3,2)));
    realBase=norm(diff(realC(:,useCams),[],2));

    cEOeuler=repmat(true,7,size(EOeuler,2));
    cEOeuler(:,baseCam)=false;
    cEOeuler(end,:)=false;
    % Lock the second camera coordinate that is largest.
    [dummy,i]=max(abs(EOeuler(1:3,2)));
    cEOeuler(i,2)=0;

    x0euler=[EOeuler(cEOeuler);OP(cOP(:))];

    ptCams=ones(1,size(markPos,2));
    cIO=logical(zeros(size(IO)));

    paramsEuler={markPos,ptCams,IO,nK,nP,EOeuler,1,OP,vis,cIO,cEOeuler,cOP};

    %z1=markPos(:,colPos(:,1));
    %z2=markPos(:,colPos(:,2));
    %[k1,a1]=convhull(z1(1,:),z1(2,:));
    %[k2,a2]=convhull(z2(1,:),z2(2,:));
    %a1=a1/prod(imSz);
    %a2=a2/prod(imSz);
    %disp(sprintf('a1=%d%%, a2=%d%%',round([a1,a2]*100)));
    a1=nan;
    a2=nan;

    maxIter=50;
    tol=1e-3;
    alphaMin=1e-3;
    c1=0.1;
    tic
    [x,code,nIter,X,alphas]=...
        gaussn_niclas_damped('pm_eulerbundle1_f',x0euler,tol,maxIter,c1,...
                             alphaMin,paramsEuler);
    time=toc;
    %code,nIter,alphas

    res.commonPts=size(OP,2);
    res.realBase=realBase;
    res.calcBase=sqrt(sum([x(1:2);referenceBase].^2));

    res.lastStep=nan;

    if (nIter>0)
	res.lastStep=max(abs(diff(X(nnz(cEOeuler)+1:end,end-1:end),[],2)));
    end

    res.minAngle=min(ptAngles);
    res.maxAngle=max(ptAngles);
    res.meanAngle=mean(ptAngles);
    res.stdAngle=std(ptAngles);
    res.medAngle=median(ptAngles);

    res.minDist=min(ptDist);
    res.maxDist=max(ptDist);
    res.meanDist=mean(ptDist);
    res.stdDist=std(ptDist);
    res.medDist=median(ptDist);

    res.x=x;
    res.x0=x0euler;
    res.code=code;
    res.iter=nIter;
    res.alphas=alphas;
    res.X=X;

    EO=EOeuler;
    EO(cEOeuler)=x(1:nnz(cEOeuler));

    res.C2=EO(1:3,2);

    f=pm_eulerbundle1_f(x,paramsEuler{:});
    m=length(f);
    n=length(x);
    res.sigma0=sqrt((f'*f)/(m-n));

    res.a1=a1;
    res.a2=a2;

    res.time=time;

    res.type='euler';

    results=res;

    camSize=0.5*max(diff(range1([OP';EO(1:3,:)'])));

    disp('Running undamped Gauss-Newton')
    tic
    [x,code,nIter,X,alphas]=...
        gaussn_niclas_undamped('pm_eulerbundle1_f',x0euler,tol,maxIter,paramsEuler);
    time=toc
    code,nIter,alphas

    if (nIter>0)
	res.lastStep=max(abs(diff(X(nnz(cEOeuler)+1:end,end-1:end),[],2)));
    end
    
    res.x=x;
    res.code=code;
    res.iter=nIter;
    res.alphas=alphas;
    res.X=X;

    EO=EOeuler;
    EO(cEOeuler)=x(1:nnz(cEOeuler));

    res.C2=EO(1:3,2);

    f=pm_eulerbundle1_f(x,paramsEuler{:});
    m=length(f);
    n=length(x);
    res.sigma0=sqrt((f'*f)/(m-n));

    res.a1=a1;
    res.a2=a2;

    res.time=time;
    
    results(2)=res;

    disp('Running damped Gauss-Newton')
    maxIter=50
    c1=-10
    tic
    [x,code,nIter,X,alphas]=...
        gaussn('pm_eulerbundle1_f',x0euler,tol,maxIter,c1,...
                             alphaMin,paramsEuler);
    time=toc
    code,nIter,alphas

    if (nIter>0)
	res.lastStep=max(abs(diff(X(nnz(cEOeuler)+1:end,end-1:end),[],2)));
    end
    
    res.x=x;
    res.code=code;
    res.iter=nIter;
    res.alphas=alphas;
    res.X=X;

    EO=EOeuler;
    EO(cEOeuler)=x(1:nnz(cEOeuler));

    res.C2=EO(1:3,2);

    f=pm_eulerbundle1_f(x,paramsEuler{:});
    m=length(f);
    n=length(x);
    res.sigma0=sqrt((f'*f)/(m-n));

    res.a1=a1;
    res.a2=a2;

    res.time=time;

    results(3)=res;

end

results(1)=results(3);

%min(30,[cat(1,results.code),cat(1,results.iter),cat(1,results.time),cat(1,results.sigma0),cat(1,results.lastStep)./cat(1,results.calcBase)*realBase*1000])

disp('Plotting')
if (~doSim)
    if (size(OP,2)~=results(1).commonPts)
        OP=zeros(3,results(1).commonPts);
        cOP=logical(ones(size(OP)));
    end
    ptOK=logical(ones(1,results(1).commonPts));
end

if (doPlot)
    if (doSim)
        % Pre-calculate point selections for features.
        iax=cell(1,length(prob.features));
        ibx=cell(1,length(prob.features));
        for i=1:length(prob.features)
            f=prob.features{i};
            if (length(f)>1)
                [dummy,ia,ib]=intersect(f,OPid);
                iax{i}=ia;
                ibx{i}=ib;
            end
        end
    end

    methodsToPresent=1:2;
    for i=methodsToPresent
        set(figure(i),'position',[(methodsToPresent(i)-1)*600,0,560,420])
        clf;
        %set(gcf,'renderer','opengl');
        set(gca,'tag','1staxes');
    end

    movies=cell(1,4);
    EOall=cell(1,4);
    iter=cat(1,results.iter);
    maxIiter=results(3).iter+2;
    %maxIiter=max([0;iter(cat(1,results.code)==0)])+2;
    for iiter=1:maxIiter
	%cpSelect,nnz(ptOK)
	%cat(1,results.iter)'
	%cat(1,results.code)'
	figure(1)
	if (any(isnan(views(cpSelect,:))))
            [az,el]=view;
	else
            az=views(cpSelect,1);
            el=views(cpSelect,2);
 end
 for rr=methodsToPresent
     iter=min(min(iiter,results(rr).iter+1),maxIter);
     switch (results(rr).type)
       case 'euler'
         EO=EOeuler;
         cEO=cEOeuler;
         EO(cEO)=results(rr).X(1:nnz(cEO),iter);
       case 'rq'
         EO=EOrq;
         cEO=cEOrq;
         EO(cEO)=results(rr).X(1:nnz(cEO),iter);
     end
     OP(cOP)=results(rr).X(nnz(cEO)+1:end,iter);

     objId=OPid;

     T=blkdiag(pm_eulerrotmat(deg2rad(rotData(cpSelect,:))),1);
     OP=applyhomoxform(T,OP);
     % Transform cameras.
     for i=1:size(EO,2)
         % Camera center.
         CC=EO(1:3,i);
         
         switch (results(rr).type)
           case 'euler'
             % Camera orientation.
             ang=EO(4:6,i);
             
             RR=pm_eulerrotmat(ang);
           case 'rq'
             RR=reshape(EO(4:12,i),3,3);
         end

         P=RR*[eye(3),-CC];
         P(4,4)=1;
         
         % Apply transformation on camera.
         P=P*inv(T);
         
         % New camera center.
         CC=P(1:3,1:3)\(-P(1:3,4));
         
         EO(1:3,i)=CC;
         
         % New camera orientation.
         switch (results(rr).type)
           case 'euler'
             ang=derotmat3d(P(1:3,1:3));
             
             EO(4:6,i)=ang';
           case 'rq'
             RR=P(1:3,1:3);
             EO(4:12,i)=RR(:);
         end
     end
     z=EOall{rr};
     z=[z,EO];
     EOall{rr}=z;
     
     figure(rr)
     ax=findobj(gcf,'tag','1staxes');
     axes(ax)
     view(az,el);
     %clf
     if (0)
         plot3(OP(1,~isCtrl),OP(2,~isCtrl),OP(3,~isCtrl),'bx',...
               OP(1,isCtrl),OP(2,isCtrl),OP(3,isCtrl),'r^');
     else
         hold on
         h=findobj(ax,'tag','points');
         OPP=XForm*OP;
         if (iiter==1 & plotx0pts)
             plot3(OPP(1,:),OPP(2,:),OPP(3,:),'.','color',x0Col,...
                   'tag','points0')
         elseif (isempty(h))
             plot3(OPP(1,:),OPP(2,:),OPP(3,:),'b.','tag','points')
         else
             set(h,'xdata',OPP(1,:),'ydata',OPP(2,:),'zdata',OPP(3,:));
         end
         hold off
     end
     axis off
     set(gca,'pos',[0,0,1,1])
     
     nImages=length(prob.cameras);
     for ii=1:length(useCams)
         i=useCams(ii);

         % Get camera icon.
         [cam,camCol]=cameraicon(camSize*[0.11,0.08,0.04]);
         cam(:,:,3)=-cam(:,:,3);
         [m,n,p]=size(cam);
         % Make camera body gray, not blue.
         %cc=reshape(camCol,m*n,p);
         %isBody=all((cc==repmat([0,0,1],m*n,1))');
         %cc(isBody,:)=repmat([0.5,0.5,1],nnz(isBody),1);
         %camCol=reshape(cc,m,n,p);
         
         % Camera center.
         CC=EO(1:3,ii);
         
         switch (results(rr).type)
           case 'euler'
             % Camera orientation.
             ang=EO(4:6,ii);
             RR=pm_eulerrotmat(ang);
           case 'rq'
             RR=reshape(EO(4:12,ii),3,3);
         end
         
         M=RR;
         CC2=CC-camSize*0.1*det(M)*M(3,:)';
         %line([CC(1),CC2(1)],[CC(2),CC2(2)],[CC(3),CC2(3)])
         
         % Apply transformation.
         T=RR*[eye(3),-CC];
         T(4,4)=1;
         cam1=reshape(applyhomoxform(inv(T),reshape(cam,m*n,p)')',m,n,p);

         T=blkdiag(XForm,1);
         cam1=reshape(applyhomoxform(T,reshape(cam1,m*n,p)')',m,n,p);

         if (~plotCameraTrace)
             tag=sprintf('cam%d',i);
             h=findobj(gca,'type','surface','tag',tag);
             if (isempty(h))
                 hold on
                 h=surf(cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,... 
                        'tag',tag,'edgealpha',1,...
                        'facealpha',0.15);
                 hold off
             else
                 set(h,'xdata',cam1(:,:,1),'ydata',cam1(:,:,2),...
                       'zdata',cam1(:,:,3));
             end
             if (iiter>iter)
                 set(h,'facealpha',1);
             end
         else
             colormap([1,0,0;0,1,0;0.5,0.5,1])
             if (ii==1 & ~isempty(findobj(gca,'type','surface','tag','cam1')))
                 % Don't redraw first camera.
             elseif (ii==1 | iiter>iter)
                 hold on
                 h=surf(cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,... 
                        'tag',sprintf('cam%d',ii),'edgealpha',1,...
                        'facealpha',1);
                 hold off
             else
                 hold on
                 h=surf(cam1(:,:,1),cam1(:,:,2),cam1(:,:,3),camCol,... 
                        'tag',sprintf('cam%d',ii),'edgealpha',0.15,...
                        'facealpha',0.15);
                 hold off
                 z=XForm*EOall{rr}(1:3,ii:length(useCams):end);
                 if (size(z,2)>1)
                     % Last point
                     z1=z(1:3,end-1);
                     % Current point.
                     z2=z(1:3,end);
                     z=[z1,z2];
                     line(z1(1,1),z1(2,1),z1(3,1),'color','r','marker','x');
                     line(z(1,:),z(2,:),z(3,:), ...
                          'color','r','linestyle','-');
                     alpha=results(rr).alphas(iiter-1);
                     if (0)
                         % Project to show where alpha=1 would have taken
                         % camera.
                         if (alpha<1)
                             z3=z1+(z2-z1)/alpha;
                             z=[z2,z3];
                             line(z(1,end-1:end),z(2,end-1:end),z(3,end-1:end), ...
                                  'color','r','linestyle','--');
                         end
                     end
                 end
             end
         end
     end
     
     if (doSim)
         % Draw 3d lines.
         xyzAll=[];
         for i=1:length(prob.features)
             ia=iax{i};
             ib=ibx{i};
             if (length(ia)>1)
                 xyz=zeros(3,length(ia));
                 xyz(:,ia)=XForm*OP(:,ib);
                 if (~any(all(xyz==0)))
                     % Don't draw features with uncalculated points.
                     xyzAll=[xyzAll,nan*ones(3,1),xyz];
                 end
             end
         end
         h=findobj(gca,'type','line','userdata',i);
         
         if (iiter==1 & plotx0pts)
             line(xyzAll(1,:),xyzAll(2,:),xyzAll(3,:),...
                  'color',x0Col,'linestyle','-','tag','lines0');
         elseif (isempty(h))
             line(xyzAll(1,:),xyzAll(2,:),xyzAll(3,:),'userdata',i,'color','b','tag','lines')
         else
             set(h,'xdata',xyzAll(1,:),'ydata',xyzAll(2,:),'zdata',xyzAll(3,:));
         end
     end
     set(gca,'projection','perspective')
     %view(az,el);
     %view(0,0)
     axis equal
     axis vis3d
     if (any(isnan(axisSizes(cpSelect,:))))
     else
         axis(axisSizes(cpSelect,:))
     end
     %axis on
     %set(gca,'xcolor','w','ycolor','w','zcolor','w');
     
     cs=['Cameras = ',sprintf('%d ',useCams),'( ',...
         sprintf('%d ',camId(useCams)),')',sprintf('Iteration %d',iter-1)];
     %title(cs)
     
     if (use2ndaxes)
         ax=findobj(gcf,'tag','2ndaxes');
         if (isempty(ax))
             ax=axes('pos',[0,0,1,1],'visible','off','tag','2ndaxes');
         end
         axes(ax);
         s=[sprintf('Case %d (camId',cpSelect),...
            sprintf(' %d',camId(useCams)),')'];
         h=findobj(ax,'tag','probnum');
         if (~isempty(h))
             set(h,'string',s);
         else
             text(0.1,0.9,s,'horizontal','left','fontsize',16,'tag','probnum');
         end
         if (iiter==1)
             s=sprintf('Alg. %s, starting approximation.',algNames{rr});
         else
             if (iiter==results(rr).iter+1)
                 if (results(rr).code==0)
                     a=results(rr).alphas(iiter-1);
                     i=inv(a);
                     if (abs(i-round(i))<1e-4)
                         sss=strrep(rats(results(rr).alphas(iiter-1)),' ','');
                     else
                         sss=sprintf('%.2f',results(rr).alphas(iiter-1));
                     end							
                     if (algIsLineSearch(rr))
                         s= sprintf('Alg. %s. \\alpha=%s. Done after %d iterations.',...
                                    algNames{rr},sss,results(rr).iter);
                     else
                         s= sprintf('Alg. %s. \\Delta=%s. Done after %d iterations.',...
                                    algNames{rr},sss,results(rr).iter);
                     end
                 else
                     s=sprintf('Alg. %s. Fail after %d iterations.',...
                               algNames{rr},results(rr).iter);
                 end
             elseif (iiter>results(rr).iter+1)
                 if (results(rr).code==0)
                     s=sprintf('Alg. %s. Done after %d iterations.',...
                               algNames{rr},results(rr).iter);
                 else
                     s=sprintf('Alg. %s. Fail after %d iterations.',...
                               algNames{rr},results(rr).iter);
                 end
             else
                 a=results(rr).alphas(iiter-1);
                 i=inv(a);
                 if (abs(i-round(i))<1e-4)
                     sss=strrep(rats(results(rr).alphas(iiter-1)),' ','');
                 else
                     sss=sprintf('%.2f',results(rr).alphas(iiter-1));
                 end							
                 if (algIsLineSearch(rr))
                     s=sprintf('Alg. %s, iter %d, \\alpha=%s.',algNames{rr},iter-1,...
                               sss);
                 else
                     s=sprintf('Alg. %s, iter %d, \\Delta=%s.',algNames{rr},iter-1,...
                               sss);
                 end
             end
         end
         h=findobj(ax,'tag','iternum');
         if (~isempty(h))
             set(h,'string',s);
         else
             text(0.1,0.1,s,'horizontal','left','fontsize',16,'tag','iternum');
         end
     end
     if (doWrite)
         if (iter>=iiter-1)
             [iter,iiter]
             if (iiter==1)
                 clear M
                 M(:,iiter)=getframe;
             else
                 M=movies{rr};
                 M(:,iiter)=getframe;
             end
             movies{rr}=M;
         end
     end
     ax=findobj(gcf,'tag','1staxes');
     axes(ax)
 end
 if (~islogical(doWait))
     if (iiter==1)
         pause % Always pause at first iteration.
     else
         pause(doWait)
     end
 elseif (doWait)
     pause
 end
 drawnow
    end
end

if (doWrite)
    for i=1:3
        movie2avi(movies{i},sprintf('/tmp/case%dalg%d.avi',cpSelect,i),'fps',1,...
                  'Videoname',sprintf('Algorithm %d',i),'quality',100,...
                  'Compression','Indeo5');
    end
end

