function confellipse(x,y,col,ax,level,type)

% Adds a confidence ellipse of level (%) to a plot based on MC data
% confellipse(x,y,col,ax,level,type)
% type 1 gives filled area, 2 gives upper and lower bound

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<6, type=1; end

if size(x)~=size(y)
    error('x and y must be of same size')
end
[N,MC]=size(x);
phi=linspace(0,2*pi,100);
for k=1:N
    P=sig2cov([x(k,:);y(k,:)]','taumax',0);
    plot(x(k,1),y(k,1),[col,'o'],'parent',ax)
    z=[cos(phi);sin(phi)];
    r=sqrt(diag(z'*(P.R)*z)')/(-log(1-level/100));
    plot(x(k,1)+cos(phi)./r,y(k,1)+sin(phi)./r,[col,'-'],'parent',ax)
end
