function y=ncfilter(b,a,u,M)

%NCFILTER performs stable non-causal filtering
%   y=ncfilter(b,a,u)
%   A stable implementation is found by forward-backward filtering
%   techniques, where poles and zeros inside the unit circle are used
%   in the forward filter, and poles and zeros outside the unit circle
%   are used in the backward filter.
%   Here, filtfilt is used. To affect the choice of initial conditions, use
%   y=ncfilter(b,a,u,M)
%   See help filtfilt for information about M.
%
%   Examples:
%      a=poly([0.5 2]);
%      u=[zeros(10,1);1;zeros(10,1)];
%      y=ncfilter(1,a,u);
%      plot([y u])
%
%   See also: filterfilter

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if nargin<4, M=[]; end

ind=find(b~=0);
k1=b(ind(1));
z=roots(b);
ind=find(abs(z)<=1);
bf=poly(z(ind));
ind=find(abs(z)>1);
bb=poly(1./z(ind));
k2=bb(end);


z=roots(a);
ind=find(abs(z)<=1);
af=poly(z(ind));
ind=find(abs(z)>1);
ab=poly(1./z(ind));
k3=ab(end);
%af,ab,bf,bb
y=k1/k2*k3*filtfilt(bf,af,bb,ab,u,M);
