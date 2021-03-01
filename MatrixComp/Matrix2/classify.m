function [digit] = classify(q,U_r)
%CLASSIFY - Classify non-seen digit by using information from previously
%           seen training data.
%
%   Takes picture of digit q and performs classification with the
%   use of unitary matrix U_r taken from SVD of training set. Returns 
%   "digit" which is what it classifies q as.
%
%   MINIMAL WORKING EXAMPLE: Classify non-seen picture q with information
%   in unitary matrix U_r, gained from using [U_r] = train(OCRtrain,rk).
%
%   digit = classify(q,U_r);

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-09-26: Initial version .
%
% Function code starts here...

%set variable err_min to inf to be sure nothing is initially bigger than it
err_min = inf;

%loop over all digits
for i = 1:length(U_r)
    %check approximation error the unseen digit q approximated by digit i
    err = norm(q - U_r{i,1}*U_r{i,1}'*q);
    %check if smallest error so far, if true set to new minimum error
    if err < err_min
        err_min = err;
        digit = i-1;
    end
end
