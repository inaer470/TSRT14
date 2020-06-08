function contour2(x,y,z,lev,levelon,varargin)
%CONTOUR2 makes a 2D version of the contour plot
%   contour2(x,y,z,lev,levelon,varargin)
%
%   contour2 is very similar to contour, but uses 2D graphics
%   lev  vector of level curve values
%   levelon binary, nonzero means print labels

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<4;
   C=contourc(x,y,z);
else
   C=contourc(x,y,z,lev);
end
if nargin<5;
   levelon=0;
end
hon=ishold;
while ~isempty(C)
   l=C(1,1); % level
   k=C(2,1); % number
   x=C(1,2:k+1);
   y=C(2,2:k+1);
   plot(x,y,varargin{:})
   hold on
   if levelon
      text(x(round(length(x)/2)),y(round(length(y)/2)),num2str(l,'%.2g'))
   end
   C=C(:,k+2:end);
end
if ~hon; hold off; end
