function [ix,newBase]=pindex(length,base)
%PINDEX Create an index vector into a parameter vector.
%
%[ix,newBase]=pindex(length,base)
%length  - length of index vector.
%base    - first index of index vector.
%ix      - index vector base:base+length-1
%newBase - base+length, i.e. first index after this vector.

% v1.0  2003-08-29. Niclas Borlin, niclas@cs.umu.se.

ix=base:base+length-1;
newBase=base+length;
