function im=jpegdecompress(coeff,quality)
%JPEGDECOMPRESS - Takes a matrix of DCT coefficients for a compressed image
%                 and a quality number (0<quality<1). Returns a
%                 decompressed image from the DCT coefficients.
%
%   MINIMAL WORKING EXAMPLE: Decompress an image from DCT coefficients
%                            using a quality number equal to 0.8:
%
%   >> q = 0.8; %set quality
%   >> imDecomp = jpegdecompress(coeff, q); %decomp image from DCT coeffs

% Author: Gustav Nystedt , guny0007@ad.umu.se
% 2018-10-24: Initial version .
%
% Function code starts here...

% Get the size of the input image.
[r,c,b]=size(coeff);
 
% Preallocate output image.
YCbCr=zeros(r,c,b,'uint8');

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

    % Restore pixels
    for i=1:8:r
        for j=1:8:c
            % Extract tile of the quantized coefficients.  
            tile=coeff(i:i+7,j:j+7,b);
            
            % De-quantize and do the inverse DCT.
            % <your code here>
            tile = (D'*(tile.*Qq)*D);
            
            % Shift back to signed, convert to uint8 and store in
            % the correct place.
            tile=tile+128;
            YCbCr(i:i+7,j:j+7,b)=uint8(tile);
        end
    end
end

% Convert image back to rgb
im=ycbcr2rgb(YCbCr);
