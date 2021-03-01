function coeff = jpegcompress(im,quality)
%JPEGCOMPRESS - Takes an RGB image and a quality number (0<quality<1) and
%               compresses it, returning the quantized DCT coefficient for
%               further decompression usage.
%
%   MINIMAL WORKING EXAMPLE: Load and compress an image in RGB format and
%                            compress it with quality number equal to 0.8:
%
%   >> im = imread(image_filename); %read image
%   >> q = 0.8; %set quality
%   >> coeff = jpegcompress(im, q); %compress image and get the coeffs

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-24: Initial version .
%
% Function code starts here...

% Get the size of the input image.
[r,c,b]=size(im);

% Verify the input image is of the correct type.
if b~=3 || ~isa(im,'uint8')
    error('Image must be uint8 RGB');
end

% Convert from RGB to YCbCr color space.
YCbCr=rgb2ycbcr(im);

% Preallocate array for coefficients
coeff=zeros(size(YCbCr));

% Get the DCT matrix and the quantization matrices.
% <your code here>
%DCT matrix
D = dctmtx(8);
Q = {};
% Luminance- and Chrominance quantization tables.
[Q{1},Q{2}] = jpegquantmat;
Q{3} = Q{2};

% For each band...
for b=1:3
    % Get the quantization matrix for this band, properly scaled.
    % <your code here>
    Qq = Q{b}/quality;
    
    % Compute coefficients
    for i=1:8:r
        for j=1:8:c
            % Extract an 8-by-8 tile, convert to floating point and
            % shift to make signed.
            tile=double(YCbCr(i:i+7,j:j+7,b))-128;

            % Do the DCT transformation and quantization
            % <your code here>
            tile = round((D*tile*D')./Qq);

            % Store the coefficients in the correct place.
            coeff(i:i+7,j:j+7,b)=tile;
        end
    end
end


