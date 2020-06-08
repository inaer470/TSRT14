classdef betadist < pdfclass
%BETADIST defines the beta distribution beta(a,b)

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Revision: 28-Oct-2019 $


properties (SetAccess = public)
  a,b;                    % Shape and scale
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=betadist(a,b);
%Constructor
%   X=betadist(a,b)  defines the distribution with a>0, b>0
%   X=betadist       defines an empty distribution
%   Y=betadist(X)    makes a copy of the distribution with a new state

%#
%#


X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.b=[]; X.a=[];  % empty structure
elseif isa(a,'betadist') % Copy PDF with new state
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
disp(['Beta(',astr,',',bstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Beta distribution';
end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression Beta(a,b)
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
str=['Beta(',astr,',',bstr,')'];
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
mu=X.a/(X.a+X.b);
end

function mu=mean(X)
%MEAN is the expectation operator
mu=X.a/(X.a+X.b);
end


function P=var(X)
%VAR is the variance operator
a=X.a; b=X.b;
P=a*b/(a+b)^2/(a+b+1);
end

function P=cov(X)
%COV is the covariance operator
P=var(X);
end

function out=skew(X)
%SKEW is the skewness operator
a=X.a; b=X.b;
out=2*(b-a)*sqrt(a+b+1)/(a+b+2)/sqrt(a*b);
end


function out=kurt(X)
%KURT is the kurtosis operator
a=X.a; b=X.b;
out=6* (a^3-a^2*(2*b-1)+b^2*(b+1)-2*a*b*(b+2)) / ( a*b*(a+b+2)*(a+b+3) );
end



%==================
%----Estimators----
%==================

function [Y,msg]=estimate(B,X)
%ESTIMATE computes a moment based estimate of X in the beta family
%   Y=estimate(betadist,X)
%   where Z=betadist
msg=[];
m=E(X);
v=var(X);
a=m*(m*(1-m)/v-1);
b=(1-m)*(m*(1-m)/v-1);
if a<=0,
   a=1;
   msg='Warning: a<0 estimated, changed to 1.';
   if nargout<2, disp(msg), end
end
if b<=0,
   b=1;
   msg='Warning: b<0 estimated, changed to 1.';
   if nargout<2, disp(msg), end
end
Y=betadist(a,b);
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
    x=linspace(0,1,N);
end
p=zeros(size(x));
ind=find(x>=0 & x<=1);
p(ind)=x(ind).^(a-1).*(1-x(ind)).^(b-1)/beta(a,b);
p=p(:);
x=x(:);
end



%==================
%----Operators-----
%==================

% Non-default operators


end
end
