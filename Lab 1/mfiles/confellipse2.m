function confellipse2(x,P,col,ax,level,linewidth)

%CONFELLIPSE2 Adds a confidence ellipse of level (%) to a SIG plot
% confellipse2(x,P,col,ax,level,linewidth)
% x is (N,2), Px is (N,2,2)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<6, linewidth=1; end
if nargin<5, level=99; end
[N,nx]=size(x);
if nx~=2
   error('nx must be 2')
end
phi=linspace(0,2*pi,100);
z=[cos(phi);sin(phi)];
%c=erfinv(ndist(0,1),level/100);
c=sqrt(erfinv(chi2dist(2),level/100));
for k=1:N
%    plot(x(k,1),x(k,2),[col,'o'],'parent',ax)
    [U,S,V]=svd(squeeze(P(k,:,:)));
    Sroot=sqrt(S);
    Ph=U*Sroot;
    plot(x(k,1)+c*Ph(1,:)*z,x(k,2)+c*Ph(2,:)*z,[col,'-'],'parent',ax,'linewidth',linewidth)
end
