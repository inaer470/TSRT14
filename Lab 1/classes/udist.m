classdef udist < pdfclass
%UDIST defines the uniform distribution U(a,b)

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $


properties (SetAccess = public)
  a,b;                    % Location parameters
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=udist(a,b);
%Constructor
%   X=udist(a,b);  % U(a,b)
%   X=udist        % Empty structure for estimate
%   Y=udist(X);    % Same distribution with different state

%#
%#

X.MC=1000; % Default number
X.xlabel={'x'};
X.state=round(1e6*rand); %Random state
if nargin==0
    X.b=[]; X.a=[];  % empty structure
elseif isa(a,'udist') % Copy PDF with new state
    Y=a;
    X.a=Y.a;
    X.b=Y.b;
    X.a=Y.a;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
elseif nargin==2
    if ~isnumericscalar(a) | ~isnumericscalar(b) | b<=a | ...
       ~isreal(a) | ~isreal(b)
        error('a and b>a must be positive scalars')
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
disp(['U(',astr,',',bstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Uniform distribution';
end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression U(a,b)
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
str=['U(',astr,',',bstr,')'];
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
mu=(X.b+X.a)/2;
end

function mu=mean(X)
%MEAN is the expectation operator
mu=(X.b+X.a)/2;
end

function P=var(X)
%VAR is the variance operator
P=(X.b-X.a)^2/12;
end

function P=cov(X)
%COV is the covariance operator
P=(X.b-X.a)^2/12;
end

function s=std(X)
%STD is the standard deviation operator
s=sqrt((X.b-X.a)^2/12);
end

function out=skew(X)
%SKEW is the skewness operator
out=0;
end

function out=kurt(X)
%KURT is the kurtosis operator
out=-6/5;
end




%==================
%----Estimators----
%==================

function [Y,msg]=estimate(U,X)
%ESTIMATE computes a moment based estimate of X in the uniform family
%   Y=estimate(U,X)
%   where U is uniform
m=E(X);
s=std(X);
b=(2*m+sqrt(12)*s)/2;
a=(2*m-sqrt(12)*s)/2;
Y=udist(a,b);
msg=[];
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
    x=linspace(X.a,X.b,N);
end
p=zeros(size(x));
ind=find(x>X.a & x<X.b);
p(ind)=1/(X.b-X.a)*ones(size(ind));
p=p(:);
x=x(:);
end

function x=rand(X,n,m)
%RAND generates random numbers in an (n,m) matrix
%   x=rand(X,n,m)  for stochastic variables X
if nargin<3; m=1; end
if nargin<2; n=1000; end
rand('state',X.state);
x=X.a+(X.b-X.a)*rand(n,m);
end

%==================
%----Operators-----
%==================

% Non-default operators

function Y=uminus(X)
%UMINUS
Y=udist(-X.b,-X.a);
end

function Y=plus(X1,X2)
%PLUS, Y=X1+X2
if isnumeric(X1)
    Y=udist(X2.a+X1,X2.b+X1);
elseif isnumeric(X2)
    Y=udist(X1.a+X2,X1.b+X2);
else % triangular distribution not supported
    MC=max([100 X1.MC X2.MC]);
    x1=rand(X1,MC,1);
    x2=rand(X2,MC,1);
    Y=empdist(x1+x2);
end
end

function Y=mtimes(X1,X2)
%MTIMES, Y=X1*X2
if isnumeric(X1)
    if X1>0
       Y=udist(X2.a*X1,X2.b*X1);
    else
       Y=udist(X2.b*X1,X2.a*X1);
    end
elseif isnumeric(X2)
    Y=udist(X1.a*X2,X1.b*X2);
else
    MC=max([100 X1.MC X2.MC]);
    x1=rand(X1,MC,1);
    x2=rand(X2,MC,1);
    Y=empdist(x1.*x2);
end
end

end
end
