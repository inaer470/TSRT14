function [yfb,ybf,M]=filtfilt(bf,af,bb,ab,u,M);

%FILTFILT Zero-phase implementation of a filter by means of forward and backward filtering.
%   y=filtfilt(b,a,u,M);          % Zero-phase non-causal filtering
%   y=filtfilt(bf,af,bb,ab,u,M);  % General non-causal filtering
%
%   For transfer functions with poles both inside and outside the unit circle,
%   a stable implementation of a filter must be non-causal.
%   The typical use is for zero-phase filters defined as |H(z)|^2 for some
%   stable TF H(z). The implementation is based on forward-backward filtering.
%
%   If an M is provided, an attempt is done to minimize the transients and
%   remove the influence of the order of application of the forward
%   and backward filter.
%   An approximation of the optimal initial conditions is used,
%   where M denotes the number of samples in both ends that are used.
%   With M=0, a default value on M is computed from the impulse response.
%
%   If different causal and non-causal filters are desired, use
%   y=filtfilt(bf,af,bb,ab,u,M);
%   Here, bf/af is the forward filter and bb/ab is the backward filter.
%   Default is bb=bf and ab=af;
%
%   The effect of transients can be judged by comparing the forward-backward
%   and backward-forward filtered sequences
%   [yfb,ybf]=filtfilt(b,a,u,M);
%
%   Example:
%     u=[5:0.1:10, 10:-0.1:2]'
%     [b,a]=getfilter(4,0.5);
%     y=filtfilt(b,a,u);
%     plot([u y])

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if nargin==5; M=[]; end;
if nargin==4; u=bb; M=ab; ab=af; bb=bf; end;
if nargin==3; u=bb; M=[]; ab=af; bb=bf; end;

u=u(:);
[N,dum]=size(u);


dum=filter(bf,af,u);
dum=filter(bb,ab,dum(N:-1:1));
yfb=dum(N:-1:1);
if nargout>1 | ~isempty(M);
  dum=filter(bb,ab,u(N:-1:1));
  ybf=filter(bf,af,dum(N:-1:1));
end;


if ~isempty(M);   % Determine the initial states and compensate for it
  if ~isnumericscalar(M) | M<=0 | round(M)~=M
     error('M must be a positive integer')
  end
  yfb0=yfb;
  ybf0=ybf;
  af=af(:);
  bf=bf(:);
  ab=ab(:);
  bb=bb(:);
  [naf,dum]=size(af);
  [nbf,dum]=size(bf);
  [nab,dum]=size(ab);
  [nbb,dum]=size(bb);
  nb=max(nab,nbb)-1;
  nf=max(naf,nbf)-1;
  if M==0;    % Determine a good M
    hf=abs(dimpulse(bf',af'));
    ind=find(hf/max(hf)>0.1);
    Mf=ind(length(ind))+1;
    hb=abs(dimpulse(bb',ab'));
    ind=find(hb/max(hb)>0.1);
    Mb=ind(length(ind))+1;
    M=max([Mf Mb]);
    % Guarantees a small impulse response: h(k)<0.1, k>M
  end;
  M=min([M N/2]);
  Ff=[ [-af(2:naf);zeros(nf-naf,1)] [eye(nf-1);zeros(1,nf-1)] ]; % observer form used in filter
  Hf=[1 zeros(1,nf-1)];
  Gf=[bf(2:nbf);zeros(nf-nbf,1)]-bf(1)*[-af(2:naf);zeros(nf-naf,1)];
  Df=bf(1);
  Of=Hf; for i=1:M-1; Of=[Of; Of(i,:)*Ff]; end;
  hf=[Df;Of(1:M-1,:)*Gf];
  St=hf'*Of; Sfb=St;
  for t=N:-1:N-M+2; St=St*Ff;          Sfb=[St;Sfb]; end;
  if af==ab & bf==bb;
    Ob=Of;
    Sbf=Sfb;
  else;
    Fb=[ [-ab(2:nab);zeros(nb-nab,1)] [eye(nb-1);zeros(1,nb-1)] ];  % observer form used in filter
    Hb=[1 zeros(1,nb-1)];
    Gb=[bb(2:nbb);zeros(nb-nbb,1)]-bb(1)*[-ab(2:nab);zeros(nb-nab,1)];
    Db=bb(1);
    Ob=Hb; for i=1:N-1; Ob=[Ob; Ob(i,:)*Fb]; end;
    hb=[Db;Ob(1:N-1,:)*Gb];
    St=hb'*Ob; Sbf=St;
    for t=N:-1:2; St=St*Fb-hb(t)*Ob(N,:); Sbb=[St;Sbf]; end;
  end;
  e=(ybf0-yfb0);
  Sfb=[zeros(M,nf); Sfb];
  Sbf=[zeros(M,nf); Sbf];
  Of=[Of; zeros(M,nf)];
  Ob=[Ob; zeros(M,nf)];
  x0Nhat= [Sfb(2*M:-1:1,:)-Of, Ob(2*M:-1:1,:)-Sbf] \ e([1:M N-M+1:N]);

  yfb=yfb0+[Sbf(2*M:-1:M+1,:)*x0Nhat(1:nf); zeros(N-2*M,1); Ob(M:-1:1,:)*x0Nhat(nf+1:2*nf)];
  ybf=ybf0+[Of(1:M,:)*x0Nhat(1:nf); zeros(N-2*M,1); Sfb(M+1:2*M,:)*x0Nhat(nf+1:2*nf)];
end;
