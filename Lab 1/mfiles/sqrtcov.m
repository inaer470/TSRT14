function Ph=sqrtcov(P)
%SQRTCOV gives a valid square root of a covariance matrix
%   Ph=sqrtcov(P)
%
%   Ph satisfies P=Ph*Ph'
%   The user has to check whether P is a valid covariance matrix with iscov
%
%   Example:
%     Generate 100 random vectors with covariance matrix P
%     x=sqrt(Ph)*randn(size(Ph,1),100);
%   See also: sqrtm

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


[U,D,V]=svd(P);
Ph=U*diag(sqrt(diag(D)));
if trace(Ph)<0
   Ph=-Ph;  % might look nicer
end
