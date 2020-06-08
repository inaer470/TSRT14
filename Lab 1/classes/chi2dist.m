classdef chi2dist < pdfclass
%CHI2DIST defines the chi2 distribution chi2(d)

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Revision: 28-Oct-2019 $


properties (SetAccess = public)
  d;                      % Degrees of freedom
  MC=1000;                % Number of Monte Carlo simulations when needed
  state=round(1e6*rand);  % Random number generator state
  xlabel='x';             % Label
end

methods
function X=chi2dist(d)
%Constructor
%   X=chi2dist(d)    defines the distribution with d being a positive integer
%   X=chi2dist       defines an exmpty distribution
%   Y=chi2dist(X)    makes a copy of the distribution with a new state

%#
%#

if nargin==0
    X.d=[];
elseif isa(d,'chi2dist') % Copy PDF with new state
    Y=d;
    X.d=Y.d;
    X.MC=Y.MC;
    X.xlabel=Y.xlabel;
    X.state=round(1e6*rand); % New unique random state
elseif nargin==1
    if ~isnumericscalar(d) | d<=0 | round(d)~=d | ~isreal(d)
        error('d must be a positive integer')
    else
	X.d=d;
    end
    X.state=round(1e6*rand); % Set unique random state
else
    error('Syntax: chi2dist(d) or chi2dist')
end
end

function X=set.d(X,d)
  if ~isnumericscalar(d) | d<=0 | round(d)~=d | ~isreal(d)
      error('d must be a positive integer')
  end
  X.d=d;
end


function disp(X, format)
%DISP
%   disp(X,format)
%   format='%11.2g'; by default
if nargin<2
   format='%11.2g';
end
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
disp(['chi2(',dstr,')'])
end

function out=desc
%DESC Returns a description of the distribution defined by the class.
out='Chi2 distribution';
end


function str=symbolic(X,format)
%SYMBOLIC returns a symbolic expression chi2(d)
%   str=symbolic(X)
%   format='%11.2g' by default
if nargin<2, format='%11.2g'; end
if isempty(X.d)
    dstr='d';
else
    dstr=num2str(X.d,format);
end
str=['chi2(',dstr,')'];
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
mu=X.d;
end

function mu=mean(X)
%MEAN is the expectation operator
mu=X.d;
end

function P=var(X)
%VAR is the variance operator
P=2*X.d;
end
function P=cov(X)
%COV is the variance operator
P=var(X);
end
function s=std(X)
%STD is the standard deviation operator
s=sqrt(2*X.d);
end
function out=skew(X)
%SKEW is the skewness operator
out=sqrt(8/X.d);
end
function out=kurt(X)
%KURT is the kurtosis operator
out=12/X.d;
end



%==================
%----Estimators----
%==================
function [Y,msg]=estimate(Z,X)
%ESTIMATE computes a moment based estimate of X in the chi2 class
%   Y=estimate(chi2dist,X)
d=round(E(X));
msg=[];
if d<=0,
   d=1;
   msg='Warning: d<0 estimated, changed to 1.';
   if nargout<2,disp(msg), end
end
Y=chi2dist(d);
end

%==================
%----Detectors-----
%==================

function [b,level,h]=detect(Z,y,pfa)
%DETECT evaluates a hypothesis test H0: y~chi2dist(ny)
%   [b,level,h]=detect(chi2dist(ny),y,pfa)
%   y     data sample
%   pfa   false alarm rate (default 0.01)
%   b     binary decision
%   level level of the test
%   h     threshold corresponding to pfa
%
%   Example:
%   Z=chi2dist(5);
%   y=rand(Z,1)
%   [b,l]=detect(Z,y,0.01)

if nargin<3, pfa=0.01; end
if nargin<2, error('CHI2DIST.DETECT: a second data argument y required'), end

if isa(y,'sig'); y=y.y; end
if length(y)>1; error('CHI2DIST.DETECT: y must be a scalar '), end
h=erfinv(Z,1-pfa);
b=(y>h);
level=erf(Z,y);
end

function roc(Z,m,h)
%ROC plots the receiver operating characteristics curve Pd as a function of Pfa
%   [pd,pfa]=roc(Z,m,h)
%   Test the hypotheses
%      H0: y~chi2dist(ny)
%      H1: y~nchi2dist(ny,m)
%   Detection formula: Pd=1-erf(erfinv(1-Pfa)-m)
%
%   With two input arguments, h is gridded automatically
%   Without output arguments, a plot is generated
%
%   Example:
%   Z=chi2dist(5);
%   roc(Z,10,0:1:20)

if nargin<3
   N=100;
   pfa=logspace(-3,0,N+1);
   pfa(end)=[];
   for i=1:N
      h(i)=erfinv(Z,1-pfa(i));
   end
else
   N=length(h);
   for i=1:N
      pfa=1-erf(Z,h);
   end
end
Z1=ncchi2dist(Z.d,m);
for i=1:N
   pd(i)=1-erf(Z1,h(i));
end
if nargout==0
   ind=1:round(N/5):N;
   plot(pfa,pd,'-',pfa(ind),pd(ind),'.')
   hold on
   for i=1:length(ind);
      text(pfa(ind(i)),pd(ind(i)),num2str(h(ind(i)),'%.3g'))
   end
   hold off
   ylabel('Pd');
   xlabel('Pfa');
end
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
    N=10000;
    x=linspace(0,X.d+15*sqrt(X.d),N);
end
p=zeros(size(x));
ind=find(x>0);
c=1/2^(X.d/2)/gamma(X.d/2);
p(ind)=x(ind).^(X.d/2-1).*exp(-x(ind)/2).*c;
p=p(:);
x=x(:);
end

function [P,x]=cdf2(X,x)
%CDF is the cumulative density function
%   P=cdf(X,x) for given (vector) x
%   [P,x]=pdf(X) automatically chooses an interval for x
%   P is (length(x),1)
if nargin<2 | isempty(x)
    N=1000;
    x=linspace(0,X.d+8*sqrt(X.d),N);
end
P=zeros(size(x));
ind=find(x>0);
P(ind)=gammainc(x(ind)/2,X.d/2)./gamma(X.d/2);
P=P(:);
x=x(:);
end





%==================
%----Operators-----
%==================

% Non-default operators
function Y=plus(X1,X2)
%PLUS, Y=X1+X2
if isa(X1,'chi2dist') & isa(X2,'chi2dist')
    Y=chi2dist(X1.d+X2.d);
else
    MC=100; %Minimum
    if isa(X1,'pdfclass')
        x1=rand(X1,MC,1);
        MC=max([MC,X1.MC]);
    else
        x1=X1;
    end
    if isa(X2,'pdfclass')
        x2=rand(X2,MC,1);
        MC=max([MC,X2.MC]);
    else
        x2=X2;
    end
    Y=empdist(x1+x2);
end
end

end %methods
end %classdef
