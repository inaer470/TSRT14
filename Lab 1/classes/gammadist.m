classdef gammadist < pdfclass
%GAMMADIST defines the gamma distribution Gamma(a,b)
%   a is a shape parameter
%   b is a scale parameter
%   Reference: http://en.wikipedia.org/wiki/Gamma_distribution

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
  a,b;                    % Shape and scale
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=gammadist(a,b);
%Constructor
%   X=gammadist(a,b)  defines the distribution with a>0, b>0
%   X=gammadist       defines an exmpty distribution
%   Y=gammadist(X)    makes a copy of the distribution with a new state

%#
%#

X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.b=[]; X.a=[];  % empty structure
elseif isa(a,'gammadist') % Copy PDF with new state
    Y=a;
    X.a=Y.a;
    X.b=Y.b;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
elseif nargin==2
    if ~isnumericscalar(a) | ~isnumericscalar(b) | a<=0 | b<=0 | ...
       ~isreal(a) | ~isreal(b)
        error('a and b must be positive scalars')
    else
        X.a=a;
        X.b=b;
    end
else
    error('a and b must be specified')
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


function disp(X,format)
if nargin<2
    format='%11.2g';
end
if isempty(X.a)
    astr='a';
else
    astr=num2str(X.a,format);
end
if isempty(X.b)
    bstr='b';
else
    bstr=num2str(X.b,format);
end
disp(['Gamma(',astr,',',bstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Gamma distribution';
end


function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression Gamma(a,b)
%   str=symbolic(X,format)
%   format is the numeric format used in mat2str
%   Default format is '%11.2g'
if nargin<2, format='%11.2g'; end
if isempty(X.a)
    astr='a';
else
    astr=num2str(X.a,format);
end
if isempty(X.b)
    bstr='b';
else
    bstr=num2str(X.b,format);
end
str=['Gamma(',astr,',',bstr,')'];
end

function n=length(X)
%LENGTH is the length of the stochastic vector X
n=1;
end


%==================
%----Moments-------
%==================

function mu=E(X)
%E is the expectation operator
mu=X.b*X.a;
end

function mu=mean(X)
%MEAN is the expectation operator
mu=X.b*X.a;
end


function P=var(X)
%VAR is the variance operator
P=X.b^2*X.a;
end

function P=cov(X)
%COV is the covariance operator
P=var(X);
end

function out=skew(X)
%SKEW is the skewness operator
out=2/sqrt(X.a);
end

function out=kurt(X)
%KURT is the kurtosis operator
out=6/X.a;
end



%==================
%----Estimators----
%==================

function Gpost=posterior(G,lh,y)
%POSTERIOR computes posterior Gamma distribution using inference from E(x)
%   Gpost=posterior(G,lh,y)
%   G   is the prior Gamma distribution
%   lh  is the likelihood distribution p(y|x),
%       for the Gamma distribution, lh can be either expdist or gammadist
%   x   is a data sampel, a vector of samples, or a sig object
%   Symbolically, the inference is computed as
%      Gpost(x;thetapost)=p(y|x)*G(x;theta)

if ~isa(lh,'expdist')| ~isa(lh,'gammadist')
%   error('GAMMADIST.POSTERIOR: lh must be either expdist or gammadist')
end

if isa(y,'sig')
   y=y.y;
end
if isa(lh,'expdist')
   Gpost=gammadist(G.a+length(y),1/(1/G.b+sum(y)));
elseif isa(lh,'gammadist')
   Gpost=gammadist(G.a+length(y)*lh.a,1/(1/G.b+sum(y)));
end
end



function [Y,msg]=estimate(G,X)
%ESTIMATE computes a moment based estimate of X
%   Y=estimate(gammadist,X)

msg=[];
m=E(X);
v=var(X);
b=v/m;
a=m/b;
if a<=0,
   a=1;
   msg='gammadist warning: a<0 estimated, changed to 1.';
   if nargout<2, disp(msg), end
end
if b<=0,
   b=1;
   msg='gammadist warning: b<0 estimated, changed to 1.';
   if nargout<2, disp(msg), end
end
Y=gammadist(a,b);
end


%==================
%----Stats---------
%==================


function [p,x]=pdf(X,x)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically chooses an interval for x
%   p is (length(x),1)

a=X.a; b=X.b;
m=E(X); s=std(X);
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(max([0,m-6*s]),m+6*s,N);
end
p=zeros(size(x));
ind=find(x>=0);
p(ind)=1/(b*gamma(a))*( (x(ind)/b).^(a-1).*exp(-x(ind)/b) );
p=p(:);
x=x(:);
end



%==================
%----Operators-----
%==================

% Non-default operators

function Y=plus(X1,X2)
%PLUS, Y=X1+X2
%   Gamma(a1,b)+Gamma(a2,b)=Gamma(a1+a2,b)
if isa(X1,'gammadist') & isa(X2,'gammadist') & X1.b==X2.b
    Y=gammadist(X1.a+X2.a,X1.b);
else
    MC=100;
    if isa(X1,'pdfclass')
        MC=max([MC X1.MC]);
    end
    if isa(X2,'pdfclass')
        MC=max([MC X2.MC]);
    end
    if isa(X1,'pdfclass')
        x1=rand(X1,MC,1);
    else
        x1=X1;
    end
    if isa(X2,'pdfclass')
        x2=rand(X2,MC,1);
    else
        x2=X2;
    end
    Y=empdist(x1+x2);
end
end

function Y=mrdivide(X1,X2)
%MRDIVIDE, Y=X1/X2
%   Gamma(a,b)/t=Gamma(a,b/t)
Y=rdivide(X1,X2);
end

function Y=rdivide(X1,X2)
%RDIVIDE, Y=X1/X2
%   Gamma(a,b)/t=Gamma(a,b/t)
if isnumeric(X1)
    Y=gammadist(X2.a,X2.b/X1);
elseif isnumeric(X2)
    Y=gammadist(X1.a,X1.b/X2);
else
    MC=100;
    if isa(X1,'pdfclass')
        MC=max([MC X1.MC]);
    end
    if isa(X2,'pdfclass')
        MC=max([MC X2.MC]);
    end
    if isa(X1,'pdfclass')
        x1=rand(X1,MC,1);
    else
        x1=X1;
    end
    if isa(X2,'pdfclass')
        x2=rand(X2,MC,1);
    else
        x2=X2;
    end
    Y=empdist(x1./x2);
end
end

function Y=mtimes(X1,X2)
%MTIMES, Y=X1*X2
%   t*Gamma(a,b)=Gamma(a,t*b)
if isnumeric(X1)
    Y=gammadist(X2.a,X2.b*X1);
elseif isnumeric(X2)
    Y=gammadist(X1.a,X1.b*X2);
else
    MC=100;
    if isa(X1,'pdfclass')
        MC=max([MC X1.MC]);
    end
    if isa(X2,'pdfclass')
        MC=max([MC X2.MC]);
    end
    if isa(X1,'pdfclass')
        x1=rand(X1,MC,1);
    else
        x1=X1;
    end
    if isa(X2,'pdfclass')
        x2=rand(X2,MC,1);
    else
        x2=X2;
    end
    Y=empdist(x1.*x2);
end
end


end
end
