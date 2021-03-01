function x=range1(A)
% x=range1(A)
% returns a two-row matrix, each column containing smallest and the largest
% element of the corresponding column in A.

% v1.0  93-04-05. Niclas Borlin, niclas@cs.umu.se.
% v1.1  94-05-07. Now distinguishes between vectors and matrices. 
% v1.2  94-09-16. Now works on complex matrices as well.
% v1.21 94-12-31. Niclas Borlin, niclas@cs.umu.se.
% v1.22 97-07-18. Changed name to range1 to avoid collision with stat tlbx.

if (prod(size(A))==0)
	x=[nan nan];
else
	if (isreal(A)) 
		x=[min(A); max(A)];
	else
		x=range1(real(A))+1i*range1(imag(A));
	end
end
