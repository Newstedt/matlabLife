function [U_r] = train(OCRtrain,rk)
%TRAIN Performs Singular value decomposition on training data and crops U
%      at desired/defined rk.
%
%   U_r = train(OCRtrain,rk), where OCRtrain is a set of pictures of
%   hand-drawn digits. Each picture is stored in a column-vector where
%   each element represents a pixel. Train performs Singular value
%   decomposition to OCRtrain and crops the U-matrix from the (rk+1):th
%   element => U_r.
%
%   MINIMAL WORKING EXAMPLE: Find U_r for training data in A with
%   rk = 15.
%
%   [U_r] = train(A, 15);

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

%pre-allocate U_r = the part of U that should be used when r=rk
U_r = cell(length(OCRtrain),1);

%loop over all digits
for j = 1:length(OCRtrain)
    %get U in SVD for each digit
    [U_r{j},~,~] = svd(OCRtrain{1,j},'econ');
    %limit U by r = rk
    U_r{j} = U_r{j}(:,1:rk(j));

end




