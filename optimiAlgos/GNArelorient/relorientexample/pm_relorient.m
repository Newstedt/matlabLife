function [P1,P2,X,inFront]=pm_relorient(x1,x2,K1,K2,ess)
%PM_RELORIENT Perform relative orientation of two cameras.
%
%[P1,P2,X,inFront]=pm_relorient(x1,x2,K1,K2[,ess])
%x1,x2   - 3xn matrices with homogenous coordinates of corresponding points.
%K1,K2   - camera calibration matrices for each camera.
%ess     - calculate essential matrix directly instead of via the 8-pt
%          algorithm. Default: false.
%P1      - canonical camera matrix [I 0].
%P2      - second camera matrix.
%X       - 4xn matrix with homogenous coordinates of object points.
%inFront - logical n-vector with 1's for points in front of both cameras.

if (ess)
    E=essmat8(K1\x1,K2\x2);
else
    % Estimate the fundamental matrix with the normalized 8-point algorithm.
    F=fundmat8(x1,x2,true,true);

    E=K2'*F*K1;
    [U,S,V]=svd(E);
    EE=U*diag((S(1,1)+S(2,2))/2*[1,1,0])*V';
end

[P1,P2,X,inFront]=camsfrome(E,K1\x1,K1\x2);
