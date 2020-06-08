function [x,w,ind]=resampling(x,w,type)
%RESAMPLING resamples particles in the particle filter
%   [x,w,ind]=resampling(x,w,type)
%
%   x are resampled according to the relative weights w
%   ind returns the resampled indexes
%
%   type  {'simple'}
%         'systematic'
%         'residual'
%         'stratified'
%
%   See: Resampling in particle filters
%        Jeroen D. Hol, T.B. Schon and F.Gustafsson
%        IEEE Nonlinear Statistical Signal Processing Workshop, 2006


% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<3
    type='simple';
end
N=length(w);

if strncmp(type,'sys',3)
  u=([0:N-1]+rand(1))/N;
  wc=cumsum(w);
  [dum,ind1]=sort([u wc(:)']);
  ind2=find(ind1<=N);
  ind=ind2- (0:N-1);
  xnew=x(ind,:);
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim1',4)
  u = cumprod(rand(1,N).^(1./[N:-1:1]));
  u = u(N:-1:1);
  wc = cumsum(w);
  wc=wc/wc(end);
  k=1;
  for i=1:N
     while(wc(k)<u(i))
        k=k + 1;
     end
     ind(i)=k;
  end
  xnew=x(ind,:);
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim2',4)
  u = cumprod(rand(1,N).^(1./[N:-1:1]));
  u = u(N:-1:1);
  wc = cumsum(w);
  wc=wc/wc(end);
  [dum,ind1]=sort([u wc(:)']);
  ind2=find(ind1<=N);
  ind=ind2- (0:N-1);
  xnew=x(ind,:);
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim3',4)
  % Maskell's fast linear in time algorithm
  w=w/sum(w);
  n=zeros(N,1);
  C=0;
  M=N;
  i=1;
  xnew=[];
  while i <= N
      if rand(1)>(1-w(i)/(1-C))^M
          n(i)=n(i)+1;
          xnew=[xnew;x(i)];
          M=M-1;
      else
          C=C+w(i);
          i=i+1;
      end
  end
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim4',4)
  % Generate sorted list of uniform numbers in u
  % Based on a Lemma in Johnson and Kotz [1970, chapter 18]
  % See Generating sorted lists of random numbers
  % Jon Louis and James B. Saxe, 1979
  wc = cumsum(w);
  wc=wc/wc(end);
  S=0;
  xnew=[];
  for i=1:N
      S=S-log(rand(1));
      u(i)=S;
  end;
  u=u/S;
  for i=1:N
      ind=find(wc<u(i));
      wc(ind)=[];
      xnew=[xnew(:);x(ind(:))];
  end
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim5',4)
  % Generate sorted list of uniform numbers in u
  % Algorithm presented in
  % Generating sorted lists of random numbers
  % Jon Louis and James B. Saxe, 1979
  wc = cumsum(w);
  wc=wc/wc(end);
  L=0;
  xnew=[];
  for i=N:-1:1
      L=L+log(rand(1))/i;
      u(i)=exp(L);
      ind=find(wc<u(i));
      wc(ind)=[];
      xnew=[xnew(:);x(ind(:))];
  end;
  wnew=ones(1,N)./N;
elseif strncmp(type,'sim',3)
  u = rand(N,1);
  wc = cumsum(w);
  wc=wc/wc(N);
  [dum,ind1]=sort([u;wc]);
  ind2=find(ind1<=N);
  ind=ind2-(0:N-1)';
  xnew=x(ind,:);
  wnew(1:N)=1/N;
elseif strncmp(type,'str',3)
  u=([0:N-1]+(rand(1,N)))/N;
  wc=cumsum(w);
  [dum,ind1]=sort([u wc(:)']);
  ind2=find(ind1<=N);
  ind=ind2- (0:N-1);
  xnew=x(ind,:);
  wnew=ones(1,N)./N;
elseif strncmp(type,'res',3)
  wa=N*w;
  nk=floor(wa);
  Nk=cumsum(nk);
  K=Nk(end);
  ind=zeros(1,K);
  k=1;
  for i=1:K
     while(Nk(k)<i)
        k=k + 1;
     end
     ind(i)=k;
  end
  wb=wa-nk;
  wb=wb/sum(wb);
  n=N-K;
  if n>0
    xnew=[x(:,ind), feval(rs,x,wb,n)];
  else
    xnew=x(:,ind);
  end
  wnew=ones(1,N)./N;
else
  error(['NL.PF Unknown sampling scheme: ',type])
end
x=xnew;
if isequal(size(w),size(wnew))
   w=wnew;
else
   w=wnew';
end
