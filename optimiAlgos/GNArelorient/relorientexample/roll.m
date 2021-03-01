function A=roll(v,l)
% function A=roll(v,l)
% Splits vector v into columns of length l

% v1.0  93-09-14. Niclas Borlin, niclas@cs.umu.se.

if (rem(prod(size(v)),l) ~= 0)
	error('Vector length must be divisible by l');
end

A=zeros(l,prod(size(v))/l);
A(:)=v;
