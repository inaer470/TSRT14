classdef ncchi2dist < pdfclass
%NCCHI2DIST defines the non-central chi2 distribution ncchi2(d,h)

%   Reference: http://en.wikipedia.org/wiki/Chi-square_distribution

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
  d;                      % Degrees of freedom
  h;                      % Non-centrality parameter
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods

function X=ncchi2dist(d,h);
%Constructor
%   X=ncchi2dist(d,h)  defines the distribution with d being a positive integer
%                      and h a positive real number
%   X=ncchi2dist       defines an exmpty distribution
%   Y=ncchi2dist(X)    makes a copy of the distribution with a new state

%#
%#

if nargin==0
    d=[];
elseif isa(d,'ncchi2dist') % Copy PDF with new state
    Y=d;
    X.d=Y.d;
    X.h=Y.h;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
elseif isa(d,'chi2dist') % Special case
    Y=d;
    X.d=Y.d;
    X.h=Y.h;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); %Random state
    X.h=0;
elseif nargin==1
    X.h=0;
    if ~isnumericscalar(d) | d<=0 | round(d)~=d | ~isreal(d)
        error('d must be a positive integer')
    else
       X.d=d;
    end
elseif nargin==2
    if ~isnumericscalar(d) | d<=0 | round(d)~=d | ~isreal(d)
        error('d must be a positive integer')
    else
       X.d=d;
    end
    if ~isnumericscalar(h) | h<0 | ~isreal(h)
        error('h must be a positive real number')
    else
       X.h=h;
    end
else
    error('Syntax: ncchi2dist(d), ncchi2dist(d,h) or ncchi2dist')
end
end

function X=set.d(X,d)
  if ~isnumericscalar(d) | d<=0 | round(d)~=d | ~isreal(d)
      error('d must be a positive integer')
  end
  X.d=d;
end

function X=set.h(X,h)
    if ~isnumericscalar(h) | h<0 | ~isreal(h)
        error('h must be a positive real number')
    else
       X.h=h;
    end
end


function disp(X, format)
%DISP
%   disp(X,format)
%   format='%11.2g'; by default
format='%11.2g';
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
if isempty(X.h)
    hstr=',h';
else
    hstr=[',',num2str(X.h,format)];
end
disp(['ncchi2(',dstr,hstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='non-central chi2 distribution';

end

function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression for ncchi2(d,h)
%   str=symbolic(X)
if nargin<2, format='%11.2g'; end
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
if isempty(X.h)
    hstr=',h';
else
    hstr=[',',num2str(X.h,format)];
end
str=['ncchi2(',dstr,hstr,')'];

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
mu=X.d+X.h;
end

function mu=mean(X)
%MEAN is the expectation operator
mu=E(X);
end

function P=var(X)
%VAR is the variance operator
P=2*X.d+4*X.h;
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
out=(X.d+3*X.h)*(2/(X.d+2*X.h))^(3/2);
end

function out=kurt(X)
%KURT is the kurtosis operator
out=12*(X.d+4*X.h)/(X.d+2*X.h)^2;
end



%==================
%----Estimators----
%==================

function [Y,msg]=estimate(Z,X)
%ESTIMATE computes a moment based estimate of X in the ncchi2 class
%   Y=estimate(ncchi2dist,X)
m=E(X);
v=var(X);
d=round((4*m-v)/2);
msg=[];
if d<=0,
   d=1;
   msg='Warning: d<0 estimated, changed to 1.';
   if nargout<2,disp(msg), end
end
h=m-d;
if h<=0,
   h=0;
   msg=char(msg,'Warning: h<0 estimated, changed to 0.');
   if nargout<2,disp(msg), end
end
Y=ncchi2dist(d,h);
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
    x=linspace(0,E(X)+6*std(X),N);
end
p=zeros(size(x));
ind=find(x>0);
d=X.d;
h=X.h;
if h==0; h=eps; end
p(ind)=0.5*(x(ind)./h).^(d/4-0.5).*exp(-0.5*(x(ind)+h)).*besseli(d/2-1,sqrt(h*x(ind)));
p=p(:);
x=x(:);
end


%==================
%----Operators-----
%==================

% Non-default operators

end % methods
end % classdef
