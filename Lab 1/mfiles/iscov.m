function [b,r]=iscov(P)
%ISCOV tests if the argument is a symmetric positive definite matrix
%
%   [ok,r]=iscov(P)
%   r=1 P contains NaN or Inf
%   r=2 P is not square
%   r=3 P is not symmetric
%   r=4 P is not positive semi-definite

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


b=true;
r=0;
if any(isnan(P)) | any(isinf(P))
   b=false;
   r=1;
elseif size(P,1)~=size(P,2)
   b=false;
   r=2;
elseif norm(P-P')/norm(P)>1e-10
   b=false;
   r=3;
else
   [U,D,V]=svd(P);
   if any(diag(D)<-eps)
       b=false;
       r=4;
   end
end
