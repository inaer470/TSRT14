classdef expdist < pdfclass
%EXPDIST defines the exponential distribution exp(mu)

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
  mu;                     % Scale parameter
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods
function X=expdist(mu);
%Constructor
%   X=expdist(mu)  defines the distribution with mu>0
%   X=expdist      defines an exmpty distribution
%   X=expdist(X)   makes a copy of the distribution with a new state

%#
%#

X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.mu=[];
elseif isa(mu,'expdist') % Copy PDF with new state
    Y=mu;
    X.mu=Y.mu;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    state=round(1e6*rand); %Random state
elseif nargin==1
    if ~isnumericscalar(mu) | mu<=0 | ~isreal(mu)
        error('mu must be a positive number')
    else
       X.mu=mu;
    end
else
    error('Syntax: expdist(mu) or expdist')
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
if isempty(X.mu)
    mustr='mu';
else
    mustr=num2str(X.mu,format);
end
disp(['exp(',mustr,')'])
end


function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Exponential distribution';
end

function str=symbolic(X)
%SYMBOLIC returns a symbolic expression exp(mu)
%   str=symbolic(X,format)
%   format is the numeric format used in mat2str
%   Default format is '%11.2g'
if nargin<2, format='%11.2g'; end
if isempty(X.mu)
    mustr='mu';
else
    mustr=num2str(X.mu,format);
end
str=['exp(',mustr,')'];
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
mu=X.mu;
end

function mu=mean(X)
%E is the expectation operator
mu=E(X);
end

function P=var(X)
%VAR is the variance operator
P=X.mu^2;
end

function P=cov(X)
%COV is the covariance operator
P=var(X);
end

function s=std(X)
%STD is the standard deviation operator
s=X.mu;
end

function out=skew(X)
%SKEW is the skewness operator
out=2;
end

function out=kurt(X)
%KURT is the kurtosis operator
out=6;
end

%==================
%----Estimators----
%==================

function [Y,msg]=estimate(Z,X)
%ESTIMATE computes a moment based estimate of X in the chi2 class
%   Y=estimate(expdist,X)
m=E(X);
if m<=0,
   m=1;
   msg='Warning: mu<0 estimated, changed to 1.';
   if nargout<2,disp(msg), end
end
Y=expdist(m); % No check for positivity
end


%==================
%----Stats---------
%==================


function [p,x]=pdf(X,x)
%PDF is the probability density function
%   p=pdf(X,x) for given (vector) x
%   [p,x]=pdf(X) automatically chooses an interval for x
%   p is (length(x),1)
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(0,10*sqrt(X.mu),N);
end
p=zeros(size(x));
ind=find(x>0);
p(ind)=1/X.mu*exp(-x(ind)/X.mu);
p=p(:);
x=x(:);
end


function [P,x]=cdf(X,x)
%CDF is the cumulative density function
%   P=cdf(X,x) for given (vector) x
%   [P,x]=pdf(X) automatically chooses an interval for x
%   P is (length(x),1)
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(0,10*sqrt(X.mu),N);
end
P=zeros(size(x));
ind=find(x>0);
P(ind)=1-exp(-x(ind)/X.mu);
P=P(:);
x=x(:);
end

function x=rand(X,n,m)
%RAND generates random numbers in an (n,m) matrix
%   x=rand(X,n,m)  for stochastic variables X
%   x=rand(X,n)    for stochastic vectors X, where x is an (n,length(X)) matrix

% Explicit inverse of P for exp
if nargin<3; m=1; end
if nargin<2; n=1; end
rand('state',X.state);
u=rand(1,n*m);
x=-log(u)*X.mu;
x=reshape(x,n,m);
end


%==================
%----Operators-----
%==================

% Non-default operators

function Y=mtimes(X1,X2)
%MTIMES, Y=X1*X2
%   t*exp(mu)=exp(t*mu)
if isnumeric(X1)
    Y=expdist(X2.mu*X1);
elseif isnumeric(X2)
    Y=expdist(X1.mu*X2);
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
