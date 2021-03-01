function [phi] = phiFun(Sigma,nu,z)

phi = exp(-0.5*z'*(Sigma\z))/((2*pi)^length(nu)/2*sqrt(det(Sigma)));

phi = phi';
