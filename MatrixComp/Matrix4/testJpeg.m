%TESTJPEG - Script for running relevant tests on our jpegcompress and
%           jpegdecompress functions.
%
%   MINIMAL WORKING EXAMPLE:
%
%   >> testJpeg.m

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-24: Initial version .
%
% Function code starts here...

clear all; clc;
%read the image
im = imread('mandrill.png');
%plot image for reference to future compressed/decompressed image
figure
imshow(im)
title(['Initial image - no compression done'])

%set quality (q = 1 initially)
q = 1;

%compress image and get the quantized DCT-coefficients
coeff = jpegcompress(im,q);
%use the DCT-coefficients and decompress the compressed image
imMod = jpegdecompress(coeff,q);

%plot the compressed/decompressed image in new window
figure
imshow(imMod)
title(['compressed and decompressed image, q = 1'])

%pre-allocate reference variables
tileExample = cell(3,1);
zeroEl = zeros(3,1);
uniEl = zeros(3,1);

%loop over all bands for only one tile
for b = 1:3
    %Look at the first tile
    tileExample{b} = coeff(1:8,1:8,b);
    %Number of zero elements
    zeroEl(b) = length(find(tileExample{b} == 0));
    %Number of distinct non-zero elements
    uniEl(b) = length(unique(tileExample{b}))-1;
end

%% Compression/Decompression using different qualities
%set vector of qualities
qSweep = [0.7, 0.1];
%pre-allocate cell array for storing images
imSweep = cell(length(qSweep),1);
%pre-allocate cell array for storing coefficients
coeffSweep = cell(length(qSweep),1);

%loop over all q's
for i = 1:length(qSweep)
   %compress image 
   coeffSweep{i} = jpegcompress(im,qSweep(i));
   %decompress image
   imSweep{i} = jpegdecompress(coeffSweep{i},qSweep(i));
   %plot image in new window
   figure
   imshow(imSweep{i})
   title(['compressed image using q = ', num2str(qSweep(i))]);
end