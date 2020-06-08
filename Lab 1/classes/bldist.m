classdef bldist < pdfclass
%BLDIST defines the band-limited distribution bl(k,d)
%   k>0 (integer) is the order
%   d>0 is the scale parameter (quantization step interpretation)
%   pdf defines as
%     p(x)=c(k)*sinc(pi*x/2/k/d)^(2k)
%   Property: its characteristic function is bandlimited, so
%     Phi(u)=FT(p(x))=0 for |u|>pi/2/d

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Revision: 28-Oct-2019 $


properties (SetAccess = public)
  k,d;                    % Order and scale (quant) level
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=bldist(k,d);
%Constructor
%   X=bldist(k,d)    defines the distribution with d being a positive integer
%   X=bldist         defines an exmpty distribution
%   Y=bldist(X)      makes a copy of the distribution with a new state

%#
%#

X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.k=[];
    X.d=[];
elseif isa(k,'bldist') % Copy PDF with new state
    Y=k;
    X.k=Y.k;
    X.d=Y.d;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
elseif nargin==1 | nargin==2
    if ~isnumericscalar(k) | k<=0 | round(k)~=k | ~isreal(k)
        error('k must be a positive integer')
    else
        X.k=k;
    end
    if nargin==2
        if ~isnumericscalar(d) | d<=0 | ~isreal(d)
            error('d must be a positive real number')
        else
            X.d=d;
        end
    else
        X.d=1;
    end
else
    error('Syntax: bldist(k), bldist(k,d) or bldist')
end
end

function out=fieldread(arg)
out=eval(arg);
end

function out=fieldwrite(arg1,arg2)
if ~all(size(eval(arg1))==size(arg2))
    error(['Cannot change the dimensions of ',arg1])
end
eval([arg1,'=arg2;'])
end


function disp(X, format)
if nargin<2
    format='%11.2g';
end
if isempty(X.k)
    kstr='k';
else
    kstr=num2str(X.k,format);
end
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
disp(['bl(',kstr,',',dstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Band-limited distribution';
end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression chi2(d)
%   str=symbolic(X)
if nargin<2, format='%11.2g'; end
if isempty(X.k)
    kstr='k';
else
    kstr=num2str(X.k,format);
end
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
str=['bl(',kstr,',',dstr,')'];
end

function n=length(X)
%LENGTH is the length of the stochastic vector X
n=1;
end

%==================
%----Moments-------
%==================

function mu=median(X)
%MEDIAN is the median operator
mu=0;
end

function mu=E(X)
%E is the expectation operator
mu=0;
end

function mu=mean(X)
%MEAN is the expectation operator
mu=0;
end


function P=var(X)
%VAR is the variance operator
tab=stats(X);
P=tab(2,X.k)*X.d^2;
end

function P=cov(X)
%COV is the variance operator
P=var(X);
end

function s=std(X)
%STD is the standard deviation operator
s=sqrt(var(X));
end

function out=skew(X)
%SKEW is the skewness operator
out=0;
end

function out=kurt(X)
%KURT is the kurtosis operator
tab=stats(X);
out=tab(3,X.k)/tab(2,X.k).^2- 3;
end

function tab=stats(X)
tab=[0.5103 0.3751 0.3030 0.2608;...
     NaN 1.14 1.66 2.25;...
     NaN NaN 8.27 14.50];
end

%==================
%----Estimators----
%==================

function [Y,msg]=estimate(Z,X)
%ESTIMATE computes an ad-hoc estimate of X in the bl class
%   Y=estimate(bldist,X)

k=3;
d=1;
for k=2:4
   Ztmp=bldist(k,1);
   dhat(k)=sqrt(var(X)/var(Ztmp));
   Zhat=bldist(k,dhat(k));
   if k==2
      e(k)=1/kurt(X);
   else
      e(k)=(kurt(Zhat)-kurt(X))/kurt(X);
   end
end
[dum,k]=min(abs(e(2:end)));
Y=bldist(k+1,dhat(k+1));
end


%==================
%----Stats---------
%==================

function x=rand(X,n,m)
%RAND genererates random numbers
%   x=rand(X,n,m)
%   x is a column vector of size N
%   An accept-reject algorithm is used, that uses the Cauchy
%   distribution as proposal density and and M=2

if nargin<2, n=1; end
if nargin<3, m=1; end
N=n*m;
rand('state',X.state);  % Repeatable random numbers
Q=cauchydist(0,1);
S=2;
M=2;
y=[];
Xn=bldist(X.k,1);
while length(y)<N
      Nx=N-length(y);
      x=rand(Q,S*Nx);
      u=rand(S*Nx,1);
      ind=find(u<pdf(Xn,x)./(M*pdf(Q,x)));
      y=[y;x(ind)];
end
x=y(1:N)*X.d;
x=reshape(x,n,m);
end


function [p,x]=pdf(X,x)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically chooses an interval for x
%   p is (length(x),1)

k=X.k;
d=X.d;
tab=stats(X);
ck=tab(1,k);

if nargin<2 | isempty(x)
    N=10000;
    x=linspace(-10,10,N);
end
x=x(:);
p=ck*ones(size(x));
z=pi*x/2/k/d;
ind=find(z~=0);
p(ind)=ck*(sin(z(ind))./z(ind)).^(2*k);
end



%==================
%----Operators-----
%==================

% Non-default operators

end
end
