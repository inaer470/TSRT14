classdef empdist < pdfclass

%EMPDIST is the empirical distribution class
%   Used for converting data vectors to a pdfclass object

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Revision: 28-Oct-2019 $

  properties (SetAccess = private)
    x;
  end

  properties (SetAccess = public)
    MC;
    state;
    xlabel;
  end

methods (Access = public)

function X=empdist(data)
%Constructor
%   X=empdist(x)  defines the empircal distribution from data vector/matrix x
%
%  Example:
%      x=randn(100,2);
%      X=empdist(x);
%      plot2(X)


MC=1000; % Default number
xlabel='';
state=round(1e6*rand); %Random state
if nargin==1 & isa(data,'empdist') % Copy PDF with new state
    X.x=data.x;
    X.MC=data.MC;
    X.state=round(1e6*rand); %Random state
elseif nargin==1
    X.state=round(1e6*rand); %Random state
    if isa(data,'sig')
       X.xlabel=data.ylabel;  % Inherit signal names
       data=data.y;         % Pick out signal samples
    else
       for i=1:size(data,2)
          X.xlabel{i}=['x',num2str(i)];
       end
    end
    if ~isnumeric(data)
        error('Data must be numeric')
    end
    X.x=data;
    X.MC=size(data,1);
else
    error('Data must be specified')
end
end



function out=fieldread(X,arg)
out=eval(X.arg);
end

function out=fieldwrite(X,arg1,arg2)
if strcmp(arg1,'x')
   error(['Cannot change the empirical data in EMP'])
else
   eval([arg1,'=',mat2str(arg2),';'])
end
end

function Xi=Xsubsref(X,S)
%SUBSREF is used to pick out parts of stochastic vector
%   Xi=X(i)
%   i is the row index/indices
S.type
if S.type~='()'
  error('EMPDIST: Use parenthesis for indexing')
end
if length(S.subs)>1
  error('EMPDIST: object is one-dimensional')
end
i=S.subs{1};
if any(i<1) | any(i>size(X.x,2))
    error('Index out of range')
end
Xi=empdist(X.x(:,i));
for k=1:length(i)
   Xi.xlabel{k}=X.xlabel{i(k)};
end
end

function disp(X)
disp(['Empirical data vector of size ',num2str(size(X.x,2)),' with ',num2str(size(X.x,1)),' samples'])
end

function str=symbolic(X)
str=['Empirical Dist'];
end

function n=length(X)
%LENGTH is the length of the stochastic vector X
n=size(X.x,2);
end

function out=desc(X)
%DESC Returns a description of the distribution defined by the class.
out='Empirical distribution';
end

%==================
%----Moments-------
%==================

function mu=median(X)
%MEDIAN	is the median operator
mu=median(X.x,1).';
end

function mu=E(X)
%E is the expectation operator
mu=mean(X.x,1).';
end

function mu=mean(X)
%MEAN is the expectation operator
mu=E(X);
end

function P=var(X)
%VAR is the variance operator
if length(X)>1
   error('X is multi-dimensional, use cov instead')
end
P=mean(X.x.^2)-mean(X.x).^2;
end


function s=std(X)
%STD is the standard deviation operator
if length(X)>1
   error('X is multi-dimensional, use cov instead')
end
s=sqrt(mean(X.x.^2)-mean(X.x).^2);
end


function P=cov(X)
%VAR is the variance operator
N=size(X.x,1);
P=0;
m=E(X);
for k=1:N;
   P=P+(X.x(k,:)-m')'*(X.x(k,:)-m');
end
P=P+P';
P=P/(2*N);
end


function out=skew(X)
%SKEW is the skewness operator
if length(X)>1
   error('skew(X) not defined for multi-dimensional variables')
end
m=ones(size(X.x,1),1)*E(X).';
m1=mean(X.x);
m2=mean((X.x-m).^2);
m3=mean((X.x-m).^3);
out=m3/m2^(3/2);
end

function out=kurt(X)
%KURT is the kurtosis operator
if length(X)>1
   error('kurt(X) not defined for multi-dimensional variables')
end
m1=mean(X.x);
m2=mean((X.x-m1).^2);
m4=mean((X.x-m1).^4);
out=m4/m2^2-3;
end


%==================
%----Estimate------
%==================

function [dist,Yhat,d]=estimate(X,pdflist)
%ESTIMATE estimates a parametric density function from a list of PDF objects
%   dist=estimate(X,pdflist)
%   Example:
%     l=list(pdfclass);
%     X=empdist(rand(100,1));
%     dist=estimate(X,l)

if ~isa(X,'empdist')
    error('X must be an empirical distribution')
end
if isstr(pdflist)
    pdflist={pdflist};
end
n=length(pdflist);
[P,x]=cdf(X);
xsort=sort(X.x);
Pemp=(1:length(xsort))'/length(xsort);
for k=1:n;
   Y=eval(pdflist{k});
   if ~isa(Y,'pdfclass')
       error([pdflist{k},' is not a valid PDF'])
   end
   Yhat{k}=estimate(Y,X);
   Phat=cdf(Yhat{k},xsort);
   e=Pemp-Phat;
   d(k)=e'*e/length(e);
end
[dum,kmin]=min(d);
dist=Yhat{kmin};
if nargout==0
   format='%11.2g';
   disp('  Norm    Distribution')
   for k=1:n
       disp([num2str(d(k),format),'   ',symbolic(Yhat{k})])
   end
end
end

%==================
%----Stats---------
%==================

function x=rand(X,n,m)
%RAND generates random numbers in an (n,m) matrix
%   x=rand(X,n,m)  for stochastic variables X
%   x=rand(X,n)    for stochastic vectors X, where x is an (n,length(X)) matrix

if nargin<3,  m=1; end
if nargin<2,  error('No dimension specified'), end
x=X.x;
rand('state',X.state)
if length(X)==1
    N=m*n;
    ind=ceil(size(x,1)*rand(N,1));
    ind=1+mod(1:N,size(x,1));  % Deterministic
    x=reshape(x(ind),n,m);
else
    N=n;
    ind=ceil(size(x,1)*rand(N,1));
    ind=1+mod(1:N,size(x,1));  % Deterministic
    x=x(ind,:);
end
end

function I=erf(X,x)
%ERF evaluates the error function I(x)=P(X<x) numerically
%   I=erf(X,x)
%   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz

if length(X)>1
    error('ERF only defined for scalar stochastic variables')
end
[P,xg]=cdf(X);
Itmp=interp(sig(P,xg),x,'degree',1);
I=Itmp.y;
end

function x=erfinv(X,I)
%ERFINV evaluates the inverse error function I(x)=P(X<x) numerically
%   x=erfinv(X,I)
%   The function generalizes the standard erfinv function to non-Gaussian distributions.
%   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
%   That is, I(x)=P(X<x)
if length(X)>1
    error('ERF only defined for scalar stochastic variables')
end
if I>1 | I<0
    error('I must be between 0 and 1')
end
[P,xg]=cdf(X);
xtmp=interp(sig(xg,P),I,'degree',1);
x=xtmp.y;
end

function [P,x]=cdf(X,x)
%CDF is the cumulative density function
%   P=cdf(X,x) for given (vector) x
%   [P,x]=cdf(X) automatically chooses an interval for x
%   p is (length(x),length(mu))

if length(X)>1
    error('CDF only defined for scalar stochastic variables')
end
if nargin<2, x=[]; end
[N,n]=size(X.x);
[xsort,ind]=sort(X.x,1);
P=(1:N)'/N;
if isempty(x)
   x=xsort;
else
   if any(x<xsort(1)) | any(x>xsort(end))
      error('Specified x value outside observation range in CDF')
   end
   %P=interp(P,xsort,x,'degree',1);
end
end

function [p,x,s]=pdf(X,x,s)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically chooses an interval for x
%   p is (length(x),length(mu))
%   s is a smoothing parameter

if length(X)>1
    error('PDF only defined for scalar stochastic variables')
end
[N,n]=size(X.x);
xmin=min(X.x);
xmax=max(X.x);
xs=sort(X.x);
xmin=xs(ceil(0.01*N));
xmax=xs(floor(0.99*N));
if nargin<2 | isempty(x);
   x=linspace(xmin,xmax,1000)';
else
   x=x(:);
end
if nargin<3 | isempty(s)
    s=(xmax-xmin)/N*30;
end
Nx=length(x);

blocksize=Nx;
if N*Nx>0
  blocksize=min(Nx,floor(1e6/N));
end
pos=1;
p=[];
while pos<=Nx
  if pos+blocksize-1>Nx
    blocksize=Nx+1-pos;
  end
  range=pos:(pos+blocksize-1);
  p(range)=sum(1/sqrt(2*pi*s^2)*exp(-0.5/s^2*(ones(N,1)*x(range)'-X.x*ones(1,blocksize)).^2))/N;
  pos=pos+blocksize;
end
p=p.';
end

%==================
%----Operators-----
%==================

function Y=uplus(X)
%UPLUS unitary plus, Y=+X
Y=empdist(X.x);
end
function Y=uminus(X)
%UMINUS unitary minus, Y=-X
Y=empdist(-X.x);
end
function Y=plus(X1,X2)
%PLUS, Y=X1+X2
try MC=X1.MC; catch MC=X2.MC; end
if nargin<3, MC=1000; end
if isa(X1,'double')  % Y=A+X
   x2=rand(X2,MC,1);
   nx=length(X2);
   if isscalar(X1)
      x1=X1*ones(MC,length(X2));
   elseif size(X1,1)==nx & size(X1,2)==1
      x1=ones(MC,1)*X1.';
   else
      error('For X1 being numeric, it must be a column vector of same size as X2 or a scalar')
   end
elseif isa(X2,'double')  % Y=Xa
   x1=rand(X1,MC,1);
   nx=length(X1);
   if isscalar(X2)
      x2=X2*ones(MC,length(X1));
   elseif size(X2,1)==nx & size(X2,2)==1
      x2=ones(MC,1)*X2.';
   else
      error('For X2 being numeric, it must be a column vector of same size as X1 or a scalar')
   end
else
   if length(X1)~=length(X2)
      error('X1 and X2 not compatible for PLUS')
   end
   x1=rand(X1,MC,1);
   x2=rand(X2,MC,1);
end
Y=empdist(x1+x2);
end

function Y=minus(X1,X2)
%MINUS, Y=X1-X2
Y=X1+(-X2);
end

function Y=times(X1,X2)
%TIMES, elementwise multiplication Y=X1.*X2
Y=mtimes(X1,X2);
end

function Y=power(X,n)
%POWER, Y=X.^n
Y=mpower(X,n);
end

function Y=mpower(X,n)
%MPOWER, Y=X^n
try MC=X1.MC; catch MC=X2.MC; end
x=rand(X,MC,1);
Y=empdist(x.^n);
end

function Y=mtimes(X1,X2)
%MTIMES, matrix multiplication Y=X1*X2
try MC=X1.MC; catch MC=X2.MC; end
if isa(X1,'pdfclass')
    MC=max([MC,X1.MC]);
    x1=rand(X1,MC,1);
else
    x1=X1;
end
if isa(X2,'pdfclass')
    x2=rand(X2,MC,1);
else
    x2=X2;
end

if isa(X1,'double')  % Y=A*X
   if size(X1,2)>1 & size(X1,2)~=length(X2)
      error('X1 and X2 not compatible for MTIMES')
   end
   Y=empdist(X2.x*X1');
elseif isa(X2,'double') %Y=XA
   if ~isscalar(X2)
      error('X1 and X2 not compatible for MTIMES')
   end
   Y=empdist(X1.x*X2);
else  % Both emp
   if length(X1)>1 & length(X2)>1
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

function Y=mldivide(A,X)
%MLDIVIDE, solves equation system AY=X
try MC=X1.MC; catch MC=X2.MC; end
if ~isa(A,'double')  % Y=A*X
   error('First argument must be numeric to MLDIVIDE')
end
if size(A,1)~=length(X)
   error('Incorrect size of A in AY=X')
end
x=rand(X,MC,1);
for k=1:MC                              %
    y(:,k)=A\x(:,k);
esend
Y=empdist(y);
end

end
end
end
