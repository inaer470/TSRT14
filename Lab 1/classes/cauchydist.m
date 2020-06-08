classdef cauchydist < pdfclass
%CAUCHYDIST defines the cauchy distribution Cauchy(m,s)
%   m is a location parameter
%   s is a scale parameter

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
  m,s;                    % Location and scale
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=cauchydist(m,s);
%Constructor
%   X=cauchydist(m,s)  defines the distribution with m>0, s>0
%   X=cauchydist       defines an exmpty distribution
%   Y=cauchydist(X)    makes a copy of the distribution with a new state

%#
%#

X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.m=[]; X.s=[];  % empty structure
elseif isa(m,'cauchydist') % Copy PDF with new state
    Y=m;
    X.m=Y.m;
    X.s=Y.s;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
elseif nargin==2
    if ~isnumericscalar(m) | ~isnumericscalar(s) | s<=0 | ...
       ~isreal(m) | ~isreal(s)
        error('m and s must be real scalars with s>0')
    else
        X.m=m;
        X.s=s;
    end
else
    error('m and s must be specified')
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
if isempty(X.m)
    mstr='m';
else
    mstr=num2str(X.m,format);
end
if isempty(X.s)
    sstr='s';
else
    sstr=num2str(X.s,format);
end
disp(['Cauchy(',mstr,',',sstr,')'])
end

function out=desc(X)
%DESC Returns a description of the distribution defined by the class.
out='Cauchy distribution';
end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression Cauchy(m,s)
%   str=symbolic(X,format)
%   format is the numeric format used in mat2str
%   Default format is '%11.2g'
if nargin<2, format='%11.2g'; end
if isempty(X.m)
    mstr='m';
else
    mstr=num2str(X.m,format);
end
if isempty(X.s)
    sstr='s';
else
    sstr=num2str(X.s,format);
end
str=['Cauchy(',mstr,',',sstr,')'];
end

function n=length(X)
%LENGTH is the length of the stochastic vector X
n=1;
end


%==================
%----Moments-------
%==================

function mu=median(X)
%median is the median operator
mu=X.m;
end

function mu=E(X)
%E is the expectation operator
mu=Inf; % By definition
end

function mu=mean(X)
%MEAN is the expectation operator
mu=E(X);
end


function P=var(X)
%VAR is the variance operator
P=Inf;
end

function P=cov(X)
%COV is the covariance operator
P=var(X);
end

function out=skew(X)
%SKEW is the skewness operator
out=Inf;
end

function out=kurt(X)
%KURT is the kurtosis operator
out=Inf;
end



%==================
%----Estimators----
%==================

function [Y,msg]=estimate(G,X)
%ESTIMATE computes a median/cdf based estimate of X
%   Y=estimate(cauchydist,X)
%   The location parameter is computed as the median of X
%   The scale parameter is somewhat arbitrarily fixed by
%   letting erfinv(X,0.9)=erfinv(G,0.9)

msg=[];
try
   m=median(X);
catch
   m=mean(X);
end
s=1;
Y=cauchydist(m,s)
p1=erfinv(Y,0.9)
p2=erfinv(X,0.9)
s=p2/p1;
Y=cauchydist(m,s)
end


%==================
%----Stats---------
%==================

  function x=rand(X,n,m)
%RAND generates random numbers from the Cauchy distribution
  %   x=rand(X,n,m)
%   x is a column vector of N random numbers
if nargin<2; n=1; end
if nargin<3, m=1; end
N=n*m;
rand('state',X.state);  % Repeatable random numbers
u=rand(N,1);
x=X.m+X.s*tan(pi*(u-0.5));
x=reshape(x,n,m);
end


function [p,x]=pdf(X,x)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically chooses an interval for x
%   p is (length(x),1)
m=X.m; s=X.s;
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(m-10*s,m+10*s,N);
end
x=x(:);
p=s/pi./((x-m).^2+s^2);
end

function [p,x]=cdf(X,x)
%CDF is the cumulative probability density function
%   p=cdf(X,x) for given (vector) x
%   [p,x]=cdf(X) automatically chooses an interval for x
%   p is (length(x),1)
m=X.m; s=X.s;
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(m-10*s,m+10*s,N);
end
x=x(:);
p=0.5+1/pi*atan((x-m)./s);
end



function I=erf(X,x)
%ERF evaluates the error function I(x)=P(X<x) numerically
%   I=erf(X,x)
%   The error function is defined as I(x)=int_{-Inf}^{x} p(z) dz
I=0.5+1/pi*atan((x-X.m)./X.s);
end

function x=erfinv(X,I)
%ERFINV evaluates the inverse error function I(x)=P(X<x) numerically
%   x=erfinv(X,I)
%   The error function is defined as I(z)=int_{-Inf}^{x} p(z) dz
%   That is, I(x)=P(X<x)

if I>1 | I<0
   error('cauchydist.erf: 0<=I<=1')
end
x=X.m+X.s*tan(pi*(I-0.5));
end

%==================
%----Operators-----
%==================

% Non-default operators
function Y=mtimes(varargin)
%MTIMES, changes scale parameter, same as times
Y=times(varargin{:});
end

function Y=times(X1,X2,MC)
%TIMES, changes scale parameter
if nargin<3, MC=1000; end
if isa(X2,'sig')  % Exception for calls like n.*y
   Y=X2.*X1;
   return
end
if isa(X1,'double')  % Y=A.*X
   if ~isscalar(X1) | X1<=0;
      error('X1 and X2 not compatible for TIMES')
   end
   Y=cauchydist(X2.m*X1,X2.s*X1);
elseif isa(X2,'double') %Y=XA
   if ~isscalar(X2) | X2<=0;
      error('X1 and X2 not compatible for TIMES')
   end
   Y=cauchydist(X1.m*X2,X1.s*X2);
else  % Both Gaussian
   x1=rand(X1,MC);
   x2=rand(X2,MC);
   Y=empdist(x1.*x2);
end
end


function Y=uminus(X)
%UMINUS,
Y=cauchydist(-X.m,X.s);
end

function Y=minus(varargin)
%MINUS, changes location parameter
Y=plus(varargin{1},-varargin{2},varargin{3:end});
end

function Y=plus(X1,X2,MC)
%TIMES, changes location parameter
if nargin<3, MC=1000; end
if isa(X2,'sig')  % Exception for calls like n.*y
   Y=X2+X1;
   return
end
if isa(X1,'double')  % Y=A+X
   if ~isscalar(X1);
      error('X1 and X2 not compatible for PLUS')
   end
   Y=cauchydist(X2.m+X1,X2.s);
elseif isa(X2,'double') %Y=XA
   if ~isscalar(X2);
      error('X1 and X2 not compatible for PLUS')
   end
   Y=cauchydist(X1.m+X2,X1.s);
else  % Both Gaussian
   x1=rand(X1,MC);
   x2=rand(X2,MC);
   Y=empdist(x1+x2);
end
end

end
end
