function [rk] = find_rk(OCRtrain)
%FIND_RK - Finds smallest r = rk that will still make the spectral norm
%          of OCRtrain-U_r*S_r*V_r' <= 0.12, where where OCRtrain is a set 
%          of pictures of hand-drawn digits. OCRtrain = USV' and U_r, S_r 
%          and V_r are cropped versions of U,S and V.
%
%   MINIMAL WORKING EXAMPLE: Find smallest possible r =rk that will still
%   make the spectral norm <= 0.12 for training set A.
%
%   rk = find_rk(A);

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

%pre-allocate rk as zero matrix
rk = zeros(length(OCRtrain),1);
%loop over all digits
for j = 1:length(OCRtrain)
    %perform SVD on current digit matrix
    [U,S,V] = svd(OCRtrain{1,j},'econ');
    %store whole A for digit j as USV', for later comparison
    A = U*S*V';
    %store norm of whole A for digit j, for later comparison
    norm_A = norm(A);
    %reset c
    c = 1;
    %reset i
    i = 20;
    %loop while spectral norm is <= 0.12
    while c == 1 && i >= 1
        %decrease number of elements included in U,S and V
        U = U(:,1:i);
        S = S(1:i,1:i);
        V = V(:,1:i);
        %compute spectral norm
        comp = norm(A-U*S*V')/norm_A;
        %check criteria for spectral norm
        if comp > 0.12
            %set rk equal to previous i -> we should take the last i
            %resulting in comp<=0.12
            rk(j) = i+1;
            c = 0;
        end
        i = i-1;
    end
end



