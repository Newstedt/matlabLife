function [QL,QC]=jpegquantmat
%JPEGQUANTMAT JPEG quantization matrices.
%
%   [QL,QC]=JPEGQUANTMAT returns the 8-by-8 quantization matrices
%   for the luminance and chroma, respectively, that are used by
%   the JPEG image compression algorithm.
%
%   Reference: ISO/IEC 10918-1: 1993(E), annex K.1.
%              https://www.w3.org/Graphics/JPEG/itu-t81.pdf

% Author: Niclas Borlin, niclas.borlin@cs.umu.se
% First version 2017-10-03.

% Tables below copied from std_luminance_quant_tbl and
% std_chrominance_quant_tbl near line 68 in jcparam.c of
% jpegsrc.v9b.tar.gz downloaded from http://www.ijg.org/files.

% Luminance quantization table.
QL = [ 16,  11,  10,  16,  24,  40,  51,  61
       12,  12,  14,  19,  26,  58,  60,  55
       14,  13,  16,  24,  40,  57,  69,  56
       14,  17,  22,  29,  51,  87,  80,  62
       18,  22,  37,  56,  68, 109, 103,  77
       24,  35,  55,  64,  81, 104, 113,  92
       49,  64,  78,  87, 103, 121, 120, 101
       72,  92,  95,  98, 112, 100, 103,  99];

% Chrominance quantization table.
QC = [ 17,  18,  24,  47,  99,  99,  99,  99
       18,  21,  26,  66,  99,  99,  99,  99
       24,  26,  56,  99,  99,  99,  99,  99
       47,  66,  99,  99,  99,  99,  99,  99
       99,  99,  99,  99,  99,  99,  99,  99
       99,  99,  99,  99,  99,  99,  99,  99
       99,  99,  99,  99,  99,  99,  99,  99
       99,  99,  99,  99,  99,  99,  99,  99];