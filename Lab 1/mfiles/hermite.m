function [p,x,w]=hermite(n)
% HERMITE polynomials as used in the quadrature rule
% Hn(x)=(-1)^n*exp(x^2)*Dn(exp(-x^2);
% p is the vector with polynomial coefficients
% x are the roots of p, computed in a numerical robust way
% w are the corresponding weights for the quadrature formula

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

p=1;
for m=2:n+1
    p=conv(p,[2 0])-[m:-1:1].*[0 0 p(1:end-1)];
end
if nargout>1
J=zeros(n,n);
   for i=1:n-1
      J(i,i+1)=sqrt(i/2);
   end
   J=J+J';%'
   [U,D]=eig(J);
   x=sqrt(2)*diag(D);
   w=U(1,:)'.^2;
end
