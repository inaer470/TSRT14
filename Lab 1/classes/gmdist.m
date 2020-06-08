classdef gmdist < pdfclass
%GMDIST defines the Gaussian Mixture distribution sum_{i=1}^k a_i N(mu_i,P_i)
%
%   help gmdist.gmdist      gives help for the constructor
%   methods gmdist          lists all methods for this class


%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $


  properties (SetAccess = public)
    alpha;
    k;
    mu;
    P;
  end

  properties (SetAccess = public)
    MC;
    state;
    xlabel;
  end

methods

function X=gmdist(varargin);
% Constructor
%   X=gmdist(alpha,mu,P)  defines the GM distribution
%   X=gmdist(k)           defines an empty mixture with k components
%   X=gmdist(X)     makes a copy of the distribution with a new state
%   For multivariate X, component j has marginal distribution
%   sum_{i=1}^k alpha(i)*N( mu(j,i), P(j,j,i) )
%   That is, the last dimension of mu and P represents the mixture,
%   and the Gaussian distribution is a special case with k=1
%
%   Example:
%   Univariate distribution
%   X=gmdist([0.2 0.8],[0 3],[1 1])
%   x=rand(X,1000)
%   plot(X,empdist(x))
%   Multivariate distribution
%   P(:,:,1)=[2 1;1 3];
%   P(:,:,2)=[3 1;1 2];
%   X=gmdist([0.6 0.4],[1 2;5 1],P)
%   plot2(X,[-5:0.5:10],[-5:0.5:10])

%#
%#


MC=1000; % Default number
alpha=[]; mu=[]; P=[]; k=[];
xlabel='';
state=round(1e6*rand); %Random state
if nargin==0
   % Do nothing
elseif nargin==1
   if isa(varargin{1},'gmdist') % Copy PDF with new state
        X=varargin{1};
        k=X.k;
        alpha=X.alpha;
        mu=X.mu;
        P=X.P;
        MC=X.MC;
        state=round(1e6*rand); %Random state
   else
        k=varargin{1};
        if isnumericscalar(k) & k>0 & round(k)==k
            %ok
        else
            error('GMDIST usage: gmdist(k) requires k to be a positive constant')
        end
    end
elseif nargin==3
   alpha=varargin{1};
   mu=varargin{2};
   P=varargin{3};
   alpha=alpha(:);
   k=size(alpha,1);
   if ~isvector(alpha)
      error('GMDIST: alpha must be a vector')
   elseif abs(sum(alpha)-1)>1e-15
     sum(alpha)-1
      error('GMDIST: alpha must be a vector that sums to one')
   end
   if size(mu,2)~=k | ~isreal(mu)
      error('GMDIST: mu must be a real matrix with size(mu,2)==length(alpha)')
   end
   n=size(mu,1);
   if ndims(P)==2
      if size(P,2)~=k | ~isreal(P)
         error('GMDIST: P must be a real matrix with third dimension equal to length(alpha)')
      else
         PP(1,1,:)=P; P=PP;
      end
   elseif  ~isreal(P)
      error('GMDIST: P must be a real matrix')
   elseif size(P,3)~=k
      error('GMDIST: last dimension of P must be equal to length(alpha)')
   end
   if size(P,1)~=n | size(P,2)~=n
      error('GMDIST: first two dimensions of P must be equal to size(mu,1)')
   end
   for m=1:k
      if ~iscov(P(:,:,m))
         error('GMDIST: P(:,:,m) must be a symmetric positive definite matrix for each m')
      end
   end
else
    error('GMDIST usage: gmdist(alpha,mu,P), gmdist(k) or gmdist')
end
X.alpha=alpha;
X.mu=mu;
X.P=P;
X.k=k;
X.xlabel=xlabel;
X.MC=MC;
X.state=state;
end


function Xi=Xsubsref(X,S)
%function Xi=subsref(X,S)
%SUBSREF is used to pick out parts of stochastic vector
%   Xi=X(i)
%   i is the row index/indices
      if S.type~='()'
        error('NDIST: Use parenthesis for indexing')
      end
      if length(S.subs)>1
         error('NDIST: object is one-dimensional')
      end
      i=S.subs{1};
      if any(i<1) | any(i>length(X.mu))
        error('Index out of range')
      end
      Xi=ndist(X.mu(i),X.P(i,i));
      for k=1:length(i)
         Xi.xlabel{k}=X.xlabel{i(k)};
      end
end

function disp(X)
str=symbolic(X);
disp(['   ',str])
end

function out=desc(X)
%DESC Returns a description of the distribution defined by the class.
out='Gaussian Mixture distribution';
end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression N(mu,P)
%   str=symbolic(X,format)
%   format is the numeric format used in mat2str
%   Default format is '%11.2g'
if nargin<2, format='%11.2g'; end
if isempty(X.mu)
    str='sum_{i=1}^k a_i N(mu_i,P_i)';
else
    str=[];
    k=X.k;
    for m=1:k
        alphastr=mat2strformat(X.alpha(m),format);
        mustr=mat2strformat(X.mu(:,m),format);
        Pstr=mat2strformat(X.P(:,:,m),format);
        str=[str,alphastr,'N(',mustr,',',Pstr,')'];
        if m<k
          str=[str,' + '];
        end
    end
end
end

function n=length(X)
%LENGTH is the length of the stochastic vector X
n=size(X.mu,1);
end


%==================
%----Moments-------
%==================

function mu=E(X)
%E is the expectation operator
mu=mean(X);
end

function mu=mean(X)
%MEAN is the expectation operator
n=size(X.mu,1);
a(1,:)=X.alpha;
w=repmat(a,n);
mu=sum( w .* X.mu, 2);
end

function P=var(X)
%VAR is the variance operator
if size(X,1)>1
   error('X is multi-dimensional, use cov instead')
end
P=cov(X);
end

function s=std(X)
%STD is the standard deviation operator
if size(X,1)>1
   error('X is multi-dimensional, use cov instead')
end
s=sqrt(cov(X));
end

function out=skew(X)
%SKEW is the skewness operator
if size(X,1)>1
   error('skew(X) not defined for multi-dimensional variables')
end
out=0;
end

function out=kurt(X)
%KURT is the kurtosis operator
if size(X,1)>1
   error('kurt(X) not defined for multi-dimensional variables')
end
out=0;
end

function P=cov(X)
%COV is the covariance operator
n=size(X.mu,1);
a(1,1,:)=X.alpha;
w=repmat(a,n,n);
mu=mean(X);
P=sum( w .* X.P,3);
for m=1:length(X.alpha);
    P=P+X.alpha(m)*(X.mu(:,m)-mu)*(X.mu(:,m)-mu)';  % Spread of the mean term
end
end

function [xhat,jhat]=map(X,type)
%MAP is the maximum a posteriori estimate
%   type='x'     MAP of X
%   type='mode'  MAP of mode
%   type='joine' joint MAP of mode and state

if nargin<2, type='x'; end
if strncmpi(type,'x',1)
   for i=1:X.k
       xx(:,i)=fminsearch(@(x) -pdf(X,x),X.mu(:,i));
       pp(i)=pdf(X,xx(:,i));
   end
   [pmax,imax]=max(pp);
   xhat=xx(:,imax);
   jhat=[];
elseif strncmpi(type,'j',1)
   [dum,jhat]=max(X.alpha);
   xhat=X.mu(:,jhat);
elseif strncmpi(type,'m',1)
   for i=1:X.k
       pp(i)=pdf(X,X.mu(:,i));
   end
   [pmax,jhat]=max(pp);
   xhat=[];
else
   error(['gmdist.map: unrecognized mode option: ',type])
end
end

%==================
%----Estimators----
%==================

function [Y,msg]=estimate(Xgmdist,Xempdist)
%ESTIMATE adapts a Gaussian mixture to empirical data
%   Y=estimate(gmdist,X)
%
%   Z denotes the Gaussian mixture distribution
%   X contains data in an empirical distribution class empdist
%
x=Xempdist.x;
N=length(x);
alpha=Xgmdist.alpha;
mu=Xgmdist.mu;
P=Xgmdist.P;
C=length(alpha);
for i=1:C
    X{i}=ndist(mu(i),P(i));
end
for k=1:10; % No termination condition yet
    % M-step: assign mode to each sample
    for i=1:C;
        w(:,i)=alpha(i)*pdf(X{i},x);
    end
    % Alt 1: MAP estimate
    [dum,imax]=max(w,[],2);
    % Alt 2: importance sampling ?!
    wc=cumsum(w,2);
    wc=wc./ ( wc(:,end)*ones(1,C) );
    u=rand(N,1);
    for j=1:N;
        dum=find(u(j)<wc(j,:));
        imax(j)=dum(1);
    end
    % E-step: estimate GM parameters
    for i=1:C
        ind=find(imax==i);
        alpha(i)=length(ind)/N;
        mu(i)=mean(x(ind));
        P(i)=var(x(ind));
        X{i}=ndist(mu(i),P(i));
    end
end
Y=gmdist(alpha,mu,P);

%disp('gmdist.estimate not available yet')
end



%==================
%----Stats---------
%==================

function x=rand(X,n,m)
%RAND generates random numbers in an (n,m) matrix
%   x=rand(X,n,m)  for stochastic variables X
%   x=rand(X,n)    for stochastic vectors X, where x is an (n,length(X)) matrix
randn('state',X.state);
if nargin<3,  m=1; end
if nargin<2,  error('No dimension specified'), end
if length(X)>1 & m>1
    error('For multivariate GM distributions X, number of columns is the length of X and cannot be set')
end
A=cumsum(X.alpha);
for j=1:n*m
    u=rand;
    ind=find(u<=A);
    ind=ind(1);
    [U,D]=svd(X.P(:,:,ind));
    x(j,:)=X.mu(:,ind).'+randn(1,length(X.mu(:,ind)))*sqrt(D)*U';
end
if m>1
   x=reshape(x,n,m);
end

end


function I=erf(X,x)
%ERF evaluates the error function I(x)=P(X<x) numerically
%   I=erf(X,x)
%   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
if nargin<2
    error('Syntax: erf(X,x)')
end
if length(X)>1
    error('ERF only defined for scalar stochastic variables')
end
I=(erf((x-mean(X))/std(X)/sqrt(2))+1)/2;
%[P,xg]=cdf(X);
%I=interp(xg,P,x,'degree',1);
disp('gmdist.erf not available yet')
end

function x=erfinv(X,I)
%ERFINV evaluates the inverse error function I(x)=P(X<x) numerically
%   x=erfinv(X,I)
%   The function generalizes the standard erfinv function to non-Gaussian distributions.
%   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
%   That is, I(x)=P(X<x)
if nargin<2
    error('Syntax: erfinv(X,I)')
end
if length(X)>1
    error('ERF only defined for scalar stochastic variables')
end
if I>1 | I<0
    error('I must be between 0 and 1')
end
disp('gmdist.erfinv not available yet')
end

function [P,x]=cdf(X,x)
%CDF is the cumulative density function
%   P=cdf(X,x) for given (vector) x
%   [P,x]=cdf(X) automatically chooses an interval for x
%   p is (length(x),1)

if length(X)>1
    error('CDF only defined for scalar stochastic variables')
end
if nargin<2, x=[]; end
[p,xsort]=pdf(X,x);
[N,n]=size(xsort);
P=cumsum(p);
P=P-P(1);
P=P/P(end);
if isempty(x)
   x=xsort;
else
   if any(x<xsort(1)) | any(x>xsort(end))
      error('Specified x value outside observation range in CDF')
   end
   %P=interp(P,xsort,x,'degree',1);
end
end

function [p,x]=pdf(X,x)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically computes a grid for x with X.MC
%   p is (length(x),length(mu))
mu=mean(X);
P=cov(X);
nx=length(X);
if nx==1; % somewhat optimized code
    if nargin<2 | isempty(x)
        x=linspace(min(mu-6*P),max(mu+6*P),X.MC)';
    else
        if size(x,2)~=1
            error('x must be a column vector for univariate X')
        end
    end
    N=length(x);
    v=squeeze(X.P(1,1,:));
    vv=repmat(v(:).',N,1);
    m=squeeze(X.mu(1,:));
    mm=repmat(m(:).',N,1);
    xx=repmat(x,1,X.k);
    p= ( (1./sqrt(2*pi*vv).*exp(-0.5*(xx-mm).^2./vv) ) *X.alpha(:) );
elseif nx>1
    if nargin>1
        [N,m]=size(x);
        if nx~=m
            error('Number of columns in x and rows in mu must be the same')
        end
    elseif nx==2
        % Grid X.MC points in nx-dimensional space
        mu=mean(X);
        P
        for k=1:nx
            x(:,k)=X.mu(:,k)+linspace(-2*sqrt(P(k,k)),2*sqrt(P(k,k)),round(10000^(1/nx)))';
        end
    else
       error('gmdist.pdf: Autogrid yet only for 1D and 2D s.v. use pdf(X,x)')
    end
    if isreal(x)
        k=0.5;
    else
        k=1;
    end
    N=size(x,1);
    m=length(X.alpha);
    for l=1:m
        Pl=squeeze(X.P(:,:,l));
        Plinv=inv(Pl);
        Pldet=det(Pl);
        for k=1:N
            epsil=x(k,:)'-X.mu(:,l);
            Vl=epsil'*Plinv*epsil;
            pp(k,l)=1/sqrt(2*pi*Pldet)*exp(-k*Vl);
        end
    end
    p=pp*X.alpha(:);
end
end

%==================
%----Operators-----
%==================

function Y=vertcat(X1,X2)
%VERTCAT [X1;X2] is used to create multivariate normal distributions
disp('gmdist.vertcat not available yet')
%Y=gmdist([X1.mu;X2.mu],diag([diag(X1.P) diag(X2.P)]));
end

function Y=uplus(X)
%UPLUS unitary plus, Y=+X
Y=gmdist(X.alpha,X.mu,X.P);
end

function Y=uminus(X)
%UMINUS unitary minus, Y=-X
Y=gmdist(X.alpha,-X.mu,X.P);
end

function Y=plus(X1,X2)
%PLUS, Y=X1+X2
if isa(X2,'sig')  % Exception for calls like n+y
   Y=X2+X1;
   return
end
if isa(X1,'double')  % Y=A+X
   mu1=X1*ones(size(X2.mu));
   Y=gmdist(X2.alpha,X2.mu+mu1,X2.P);
elseif isa(X2,'double')  % Y=Xa
   mu2=X2*ones(size(X1.mu));
   Y=gmdist(X1.alpha,X1.mu+mu2,X1.P);
elseif isa(X2,'gmdist')
       mu1=X1.mu; P1=X1.P; a1=X1.alpha; k1=X1.k;
       mu2=X2.mu; P2=X2.P; a2=X2.alpha; k2=X2.k;
       mu=[]; P=[]; a=[];
       for i=1:k1;
         for j=1:k2;
            a=[a a1(i)*a2(j)];
            mu=[mu mu1(i)+mu2(j)];
            P=[P P1(i)+P2(j)];
         end
       end
       a=a/sum(a);
       Y=gmdist(a,mu,P);
elseif isa(X2,'ndist')
       mu1=X1.mu; P1=X1.P;
       mu2=repmat(X2.mu,[1,X1.k]);
       P2=repmat(X2.P,[1,1,X1.k]);
       if size(mu1,1)~=size(mu2,1)
          error('X1 and X2 not compatible for PLUS')
       end
       Y=gmdist(X1.alpha,mu1+mu2,P1+P2);
elseif isa(X2,'pdfclass')
       MC=max([100 X1.MC X2.MC]);
       x1=rand(X1,MC,1);
       x2=rand(X2,MC,1);
       Y=empdist(x1+x2);
else
   error('Incorrect argument to PLUS in GMDIST class')
end
end

function Y=minus(X1,X2)
%MINUS, Y=X1-X2
Y=X1+(-X2);
end

function Y=mrdivide(X1,X2,MC)
%MRDIVIDE, division Y=X1/X2
if isa(X2,'sig')  % Exception for calls like n./y
   Y=X2/X1;
elseif isa(X2,'double')
   Y=X1*(1/X2);
end
end

function Y=rdivide(X1,X2,MC)
%RDIVIDE, elementwise division Y=X1./X2
if isa(X2,'sig')  % Exception for calls like n./y
   Y=X2./X1;
elseif isa(X2,'double')
   Y=X1.*(1./X2);
end
end

function Y=times(X1,X2,MC)
%TIMES, elementwise multiplication Y=X1.*X2
if isa(X2,'sig')  % Exception for calls like n.*y
   Y=X2.*X1;
   return
end
if nargin<3, MC=1000; end
if isa(X1,'double')  % Y=A.*X
   if size(X1,1)~=size(X2.mu,1);
      error('X1 and X2 not compatible for TIMES')
   end
   Y=ndist(X1.*X2.mu,diag(X1)*X2.P*diag(X1));
elseif isa(X2,'double') %Y=XA
   if size(X1.mu,1)~=size(X2,1);
      error('X1 and X2 not compatible for TIMES')
   end
   Y=ndist(X1.mu.*X2,diag(X2)*X1.P*diag(X2));
else  % Both Gaussian
   x1=rand(X1,MC);
   x2=rand(X2,MC);
   Y=empdist(x1.*x2);
end
end

function Y=power(X,n,MC)
%POWER, Y=X.^n
if nargin<3; MC=1000; end
Y=mpower(X,n,MC);
end

function Y=mpower(X,n,MC)
%MPOWER, Y=X^n
if X.mu==0 & X.P==1 & n==2
   Y=chi2dist(1);
else
   if nargin<3; MC=1000; end
   x=rand(X,MC,1);
   Y=empdist(x.^n);
end
end

function Y=mtimes(X1,X2,MC)
%MTIMES, matrix multiplication Y=X1*X2
if nargin<3, MC=1000; end
if isa(X1,'double')  % Y=A*X
   if ~isscalar(X1) & size(X1,2)~=size(X2.mu,1);
      error('X1 and X2 not compatible for MTIMES')
   end
   Y=ndist(X1*X2.mu,X1*X2.P*X1');
elseif isa(X2,'double') %Y=XA
   if ~isscalar(X2)
      error('X1 and X2 not compatible for MTIMES')
   end
   Y=ndist(X1.mu.*X2,X2.^2*X1.P);
         %elseif isa(X1,'ndist') & isa(X2,'ndist')   % Both Gaussian
else  % Both distributions
   if length(X1)>1 | length(X2)>1
      error('X1 and X2 not compatible for MTIMES')
   end
   x1=rand(X1,MC,1);
   if length(X1)==1
      x1=x1*ones(1,length(X2));
   end
   x2=rand(X2,MC,1);
   if length(X2)==1
      x2=x2*ones(1,length(X1));
   end
   Y=empdist(x1.*x2);
end
end


%==================
%----Plots---------
%==================


function plot(varargin)
%PLOT illustrates one or more PDF objects
%   plot calls pdfplot
pdfplot(varargin{:},'view','pdf')
end


function cdfplot(varargin)
%CDFPLOT illustrates the CDF of X
%   plot calls pdfplot
pdfplot(varargin{:},'view','cdf')
end

function erfplot(varargin)
%CDFPLOT illustrates the error function ERF of X
%   plot calls pdfplot
pdfplot(varargin{:},'view','erf')
end


function contour(p,ind)
%CONTOUR illustrates the PDF object p in 2D
%   contour(p,ind)
%   ind is the two indices X(ind) to plot (default the first two ones)

if nargin<2, ind=[1 2]; end
nx=length(p);
if nx==1
    error('contour requires the s.v. X to be at least two-dimensional')
end
if length(ind)~=2;
    error('index vector ind must be two dimensional')
end
%[p,x]=pdf(X(ind));
%surf(x(:,1),x(:,2),histeq(p))
%mesh(x(:,1),x(:,2),histeq(p))
mu=E(p);
P=cov(p);
[U,S,V]=svd(P);
Sroot=sqrt(S);
Ph=U*Sroot;
N=100;
phi=linspace(0,2*pi,N);
plot(mu(1),mu(2),'+')
hold on
plot(mu(1)+Ph(1,:)*[cos(phi);sin(phi)],...
mu(2)+Ph(2,:)*[cos(phi);sin(phi)],'-')
plot(mu(1)+2*Ph(1,:)*[cos(phi);sin(phi)],...
     mu(2)+2*Ph(2,:)*[cos(phi);sin(phi)],'-')
plot(mu(1)+3*Ph(1,:)*[cos(phi);sin(phi)],...
     mu(2)+3*Ph(2,:)*[cos(phi);sin(phi)],'-')
hold off
axis('equal')
legend('mean','1 std','2 std','3 std')
end

function surf(p,ind)
%SURF illustrates the PDF object p in 2D using a surf plot
%   surf(p,ind)
%   ind is the two indices X(ind) to plot (default the first two ones)

if nargin<2, ind=[1 2]; end
nx=length(p);
if nx==1
    error('contour requires the s.v. X to be at least two-dimensional')
end
if length(ind)~=2;
    error('index vector ind must be two dimensional')
end
[p,x]=pdf(p(ind));
%surf(x(:,1),x(:,2),p)
surf(x(:,1),x(:,2),p)
end

end
end
