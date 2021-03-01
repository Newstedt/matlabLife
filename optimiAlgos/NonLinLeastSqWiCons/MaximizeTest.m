% Test maximization function on objective function
% Create data
CorrMatDat = [1 0.8; 0.8 1]; nu = [5 20]'; N = 100;
dim = length(nu);
[X] = myRandTsamp(N,nu,CorrMatDat);

% Use univar t-cdf on sample data to transform into uniformly dist data
u = zeros(N,dim);
for i = 1:dim
    u(:,i) = tcdf(X(:,i),nu(i));
end

% offSig = 0.6:0.1:0.9;
% nu = [3:7;18:22]';
% MLfuncVec = []; thetaMat = [];
% 
% for i = 1:length(offSig)
%     for j = 1:length(nu(:,1))
%         for k = 1:length(nu(:,2))
%             Sigma = [1 offSig(i) offSig(i) 1]';
%             theta = [[nu(j,1),nu(k,2)]';Sigma];
%             thetaMat = [thetaMat,theta];
%             MLfunc = objFunc1(u,dim,theta);
%             MLfuncVec = [MLfuncVec,MLfunc];
%         end
%     end
%     fprintf('ellapse %d \n',i)
% end
Sigma0 = [1 0.8; 0.8 1];
A = chol(Sigma0,'Lower');
A = reshape(A',[],1);

theta0 = [nu;A];

MLfunc = @(theta) -objFunc1(u,dim,theta);

C = []; b = []; Ceq = []; beq = [];
nonlcon = @constraintFun;
thetaMin = fmincon(MLfunc,theta0,C,b,Ce);
