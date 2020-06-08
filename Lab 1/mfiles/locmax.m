function [ymax,ind]=locmax(y,n)

%LOCMAX finds the local maxima in the vector y sorted in order
%   [ymax,ind]=locmax(y,n)
%   ymax  Local maxima
%   ind   Corresponding indexes
%   n     Return the n largest maxima (default all)
%
%   Examples:
%   y=exsig('sinc');
%   [ymax,ind]=locmax(y.y);
%   plot(y.t,y.y,'-',y.t(ind(1:6)),ymax(1:6),'o')

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


% Find all local maxima
ind=find( (diff([y(:)' 0])<0).*(diff([0 y(:)'])>0) >0);
% Sort them
[ymax,ind2]=sort(y(ind));
ind=ind(ind2(end:-1:1));
ymax=ymax(end:-1:1);
if nargin>1
   ymax=ymax(1:n);
   ind=ind(1:n);
end
