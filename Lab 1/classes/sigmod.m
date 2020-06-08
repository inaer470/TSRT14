classdef sigmod < nl

%SIGMOD is the model object for signal models
%   sigmod is a special case of the NL object with nx=0
%   m=sigmod(h,nn)
%   The signal model is defined as
%        y(t) = h(t,th) + e(t),  e(t)~pe
%   h is an inline object or a string with arguments t,th,
%   and nn=[ny nth] gives the dimensions of y and the parameters th.
%   Convention of dimensions:
%     y (ny,1) vector or (ny,N) matrix
%     th vector of length nth
%
%   The distributions pe is specified as pdfclass objects, or
%   as a covariance matrix of suitable dimension, in which
%   case a ndist object is created.
%
%   help sigmod.sigmod    gives help for the constructor
%   methods sigmod        lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

methods

function m=sigmod(h,nn)
%SIGMOD constructor signal models
%   m=sigmod(h,nn)
%   The signal model is defined as
%        y(t) = h(t,th) + e(t),  e(t)~pe
%   h is an inline object or a string with arguments t,th,
%   and nn=[ny nth] gives the dimensions of y and the parameters th.
%   Convention of dimensions:
%     y (ny,1) vector or (ny,N) matrix
%     th vector of length nth
%
%   The distributions pe is specified as pdfclass objects, or
%   as a covariance matrix of suitable dimension, in which
%   case a ndist object is created.
%   Other fields: name thlabel ylabel
%
%   Examples:
%   s=sigmod('sin(2*pi*th(1)*t)',[1 1]);
%   s.pe=0.01;  s.th=0.1;
%   s.thlabel='Frequency'; s.name='Sine wave in noise';
%   y=simulate(s,1:10);
%   shat=estimate(s,y)

%#
%#


if isa(h,'sigmod')
   m.f=h.f; m.h=h.h; m.nn=h.nn; m.fs=h.fs;
   m.pv=h.pv;
   m.pe=h.pe;
   m.x0=h.x0;
   m.th=h.th;
   m.px0=h.px0;
   m.xlabel=h.xlabel;
   m.ylabel=h.ylabel;
   m.ulabel=h.ulabel;
   m.name=h.name;
elseif nargin==1
   error('sigmod: two input arguments required')
else
  mnl=nl('x',h,[0 0 nn],NaN);
   m.f=mnl.f; m.h=mnl.h; m.nn=mnl.nn; m.fs=mnl.fs;
   m.pv=mnl.pv;
   m.pe=mnl.pe;
   m.x0=mnl.x0;
   m.th=mnl.th;
   m.px0=mnl.px0;
   m.xlabel=mnl.xlabel;
   m.ylabel=mnl.ylabel;
   m.ulabel=mnl.ulabel;
   m.name=mnl.name;
end
end

function z=simulate(s,t,varargin)
%SIMULATE evaluates s(t)
%   y=simulate(s,t)
%
%   s sigmod object
%   t time instant or vector (default 0)
%   y sig object
%
%
%   Property   Value       Description
%   -----------------------------------------------------------------
%   MC         {30}        Number of MC simulations when th uncertain
%

if nargin<2
   try
      z=simulate(s,0);
      return
   catch
      error('SIGMOD.simulate: time argument required, s(0) not defined')
   end
end
opt=struct('MC',30);
opt=optset(opt,varargin);

y=[];
N=length(t);
ny=s.nn(3);
try
    y=feval(s.h,t(:),[],[],s.th);
catch
    for k=1:N
        yk=feval(s.h,t(k),[],[],s.th);
        y=[y;yk(:)];
    end
end
if isa(s.pe,'pdfclass')
  y=y+rand(s.pe,N);
end

yMC=[];
xMC=[];
if ~isempty(s.P) & any(any(s.P~=0))
       nx=s.nn(1); nu=s.nn(2); ny=s.nn(3); nth=s.nn(4); N=length(t);
       yMC=zeros(opt.MC,N,ny);
       thMC=rand(ndist(zeros(nth,1),s.P(1:nth,1:nth)),opt.MC)';
       mtmp=s;
       for i=1:opt.MC;
           mtmp.th=s.th+thMC(:,i);
           ztmp=simulate(mtmp,t(:),varargin{:},'MC',0);
           yMC(i,:,:)=ztmp.y;
       end
end
%       size(y),size(t)
z=sig(y,t,[],[],yMC);
z.fs=s.fs;
z.MC=opt.MC;
z.ylabel=s.ylabel;
end


%==================
%----Estimators----
%==================


function [shat,thhat,P]=ls(s,y);
%LS computes the least squares estimate of theta
%   [shat,thhat,P]=ls(s,y);
%
%   The function computes
%   argmin V(x)= argmin sum( (y-H*th)^T*(y-H*th))
%   For nonlinear functions s.h, numgrad is used to
%   linearize around s.th, which then becomes a crucial parameter.
%
%   Example:

nth=s.nn(4);
ny=s.nn(3);
if isempty(s.pe)
   R=eye(ny);
else
   R=cov(s.pe);
end
Rinv=inv(R);
RN=zeros(nth);
fN=zeros(nth,1);
% code for ny=1 only
for k=1:length(y.t)
   H=numgrad(s.h,4,y.t(k),[],[],s.th);
   yhat=bsxfun(@plus, feval(s.h,y.t(k),[],[],s.th), mean(s.pe));
   RN=RN+H'*Rinv*H;
   fN=fN+H'*Rinv*(y.y(k)-yhat+H*s.th);
end
thhat=RN\fN;
P=inv(RN);
shat=s;
shat.th=thhat';
shat.P=P;
end

function [xhat,shat]=wls(s,y);
%WLS computes the weighted least squares estimate
%   [xhat,shat]=wls(s,y);
%
%   The function computes
%   argmin V(x)= argmin sum( (y-H*x)^T*inv(cov(s.pe))*(y-H*x))
%   For nonlinear functions s.h, numgrad is used to
%   linearize around s.x0, which then becomes a crucial parameter.
%
%   Example:
%   s=exsensor('toa',3,1);
%   y=simulate(s,1);
%   xhat=wls(s,y);
%   plot(s)
%   hold on
%   xplot2(xhat,'conf',90)
%   hold off


R=cov(s.pe);
Rinv=inv(R);
for k=1:length(y.t)
   H=numgrad(s.h,2,y.t(k),s.x0,y.u(k,:)',s.th);
   yhat=bsxfun(@plus, feval(s.h,y.t(k),s.x0,y.u(k,:)',s.th), mean(s.pe));
   PP=pinv(H'*Rinv*H);
   ff=H'*Rinv*(y.y(k,:)'-yhat+H*s.x0);
   x(k,:)=PP*ff;
   yy(k,:)=x(k,:)*H';
   PP=0.5*(PP+PP');
   Px(k,:,:)=PP;
   Py(k,:,:)=H*PP*H';
end
xhat=sig(yy,y.t,y.u,x,Py,Px);
if length(y.t)==1
   shat=s;
   shat.x0=x';
   shat.px0=squeeze(Px(end,:,:));
end
end

function [xhat,shat,res]=ml(s,y);
%ML computes the ML/NLS parameter estimate in the sigmod object s from y
%   [xhat,shat,res]=ml(s,y);
%
%   The function computes the NLS (ML for Gaussian s.pe) estimate
%   argmin V(x)= argmin sum( (y-s.h(x))^T*inv(cov(s.pe))*(y-s.h(x)))
%
%   s      Signal model, where s.x0 contains the initial value of x
%   y      Observations from the true signal model s0
%   shat   Estimated signal model, with shat.x0=xhat
%   xhat   Estimated parameter vector as sig object
%   res    The structure returned from the NLS function. Here res.TH gives
%          the parameter estimate after each iteration for instance
%
%   The function calls nl.estimate with thmask=zeros(nth,1)
%
%   Example:
%   s=exsensor('toa',3,1);
%   y=simulate(s,1);
%   xhat=ml(s,y);
%   plot(s)
%   hold on
%   xplot2(xhat,'conf',90)
%   hold off

nth=s.nn(4);
[shat,res]=estimate(s,y,'thmask',zeros(nth,1));
yhat=simulate(shat,y);
xhat=sig(yhat.y,y.t,y.u,yhat.x);
end

function x=crlb(s,y)
%CRLB computes the Cramer-Rao lower bound for the state at y.x
%   x=crlb(s,y)
%
%   s   is a sigmod object
%   y   is a sig object of length 1, where y.x (and y.u) is used,
%       if y is omitted, s.x0 is used as x
%   x   is a sig object where x.x=y.x and x.Px contains the CRLB
%   Without output arguments, a confidence ellipse is plotted
%
%   See also:
%     crlb2 plots trace, det or max(eig) over a 2D grid
%     lh1, lh2 for computing/plotting the likelihood for an observation y
%
%   Example:
%   s=exsensor('toa',3,1);
%   plot(s)
%   hold on
%   crlb(s)
%   hold off

if isnan(s.fs)
   s.fs=1;  % Cannot be NaN in nl.crlb
end
if nargin<2
    x0=s.x0';
    u=[];
else
  if length(y)>1
     error('CRLB can only be computed for one sample at the time')
  end
  x0=y.x';
  u=y.u';
end
y=sig(x0',s.fs,u',x0');
x=crlb(nl(s),y);
if nargout==0
   xplot2(x,'conf',90)
end
end

function [cx,X1,X2]=crlb2(s,y,x1,x2,ind,type);
%CRLB2 computes a scalar CRLB measure over a 2D state space grid
%   [cx,X1,X2]=crlb2(s,y,x1,x2,ind,type);
%
%   s   is a sigmod object
%   y   is a sig object of length 1specifying y.x.
%       If y is empty, x0 is taken from s.x0
%   ind  {[1 2]}  A vector with two integers indicating which
%                 states x1,x2 refer to in x
%   type  {'trace'}, 'rmse', 'det', 'max'
%                 operator for transforming P(x) to scalar c(x)
%                 rmse=sqrt(trace)
%
%   Symbolically, cx=type(crlb(s,x))
%
%   Example:
%   s=exsensor('toa',2,1);
%   plot(s);
%   hold on, crlb(s); hold off

if nargin<6
   type='trace';
end
if nargin<5
   ind=[1 2];
end
if nargin<4
   x2=linspace(0,1,100);
end
if nargin<3
   x1=linspace(0,1,100);
end
if length(ind)~=2
   error('SIGMOD.CRLB2: ind must be a vector of length 2')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SIGMOD.CRLB2: ind must be a vector with positive integers')
end
if  any(ind>s.nn(1))
   error('SIGMOD.CRLB2: The state indeces in ind cannot exceed nx=s.nn(1)')
end
if isempty(s.pe)
   error('SIGMOD.CRLB2: pe cannot be empty for crlb2')
end
if ~isstr(type)
   error('SIGMOD.CRLB2: type must be a string')
end
[X1,X2]=meshgrid(x1,x2);
X11=X1(:);
X22=X2(:);
N1=length(x1);
N2=length(x2);
N=N1*N2;


if isempty(y)
    x0=s.x0;
    u=[];
else
  if length(y)>1
     error('SIGMOD.CRLB2 can only be computed for one sample at the time')
  end
  x0=y.x';
  u=y.u';
end


X=repmat(x0',N,1);
X(:,ind)=[X1(:) X2(:)];
ss=s;
for i=1:N;
   xi=X(i,:)';
   x=crlb(ss,sig(xi',s.fs,u',xi'));
   P=squeeze(x.Px(1,:,:));
   switch lower(type)
   case 'rmse'
      cx(i)=sqrt(trace(P));
   case 'trace'
      cx(i)=trace(P);
   case 'det'
      cx(i)=det(P);
   case 'max'
      cx(i)=max(eig(P));
   otherwise
       error(['SIGMOD.CRLB2: unknown option type = ',type])
   end
end
cx=reshape(cx,N1,N2);
if nargout==0
   [cs,h]=contour(x1,x2,cx);
   clabel(cs,h)
   xlabel(gca,s.xlabel{ind(1)})
   ylabel(gca,s.xlabel{ind(2)})
end
end

function [lh,px,px0,x1]=lh1(s,y,x1,ind);
%LH1 computes the one-dimensional likelihood function over a state space grid
%   [lh,px,px0,x]=lh1(s,y,x,ind);
%
%   ind  {[1]}    An integer indicating which state to vary
%   x1            Vectors with grid points to evaluate lh over
%   y             Measurement(s) as sig object
%   lh            Likelihood function
%   px0           Prior evaluated in x1
%   px            Posterior=px0*lh
%
%   Symbolically, px=prod_{t=y.t} s.pe(y(t);s.h(t,x,u(t),s.th|x(ind)))
%
%   Example:
%   s=exsensor('toa',2,1);
%   y=simulate(s,1);
%   subplot(2,1,1), plot(s);
%   subplot(2,1,2), lh1(s,y);
%   See also the larger example in lh2

if nargin<4
   ind=[1];
end
if nargin<3 | isempty(x1)
   x1=linspace(0,1,100);
end
if length(ind)>1
   error('SIGMOD.LH1: ind must be an integer')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SIGMOD.LH1: ind must be a positive integer')
end
if ind>s.nn(1)
   error('SIGMOD.LH1: state index ind cannot exceed nx in s.nn(1)')
end
if isempty(s.pe)
   error('SIGMOD.LH1: pe cannot be empty for lh1')
end
N=length(x1);
X=repmat(s.x0,1,N);
X(ind,:)=x1(:)';
if ~isempty(s.px0)
   px0=pdf(s.px0,X');
else
   px0=ones(N,1);
end
px=px0;
for k=1:length(y.t);
   yhat=feval(s.h,y.t(k),X,y.u(k,:),s.th);
   if size(yhat,2)~=N
      yhat=zeros(s.nn(3),N);
      for n=1:N
         yhat(:,n)=feval(s.h,y.t(k),X(n,:)',y.u(k,:),s.th);
      end
   end
   epsi=repmat(y.y(k,:).',1,N)-yhat;
   lh =pdf(s.pe,epsi.');
   px=px.*lh;
end
if nargout==0
   plot(x1,lh)
   xlabel(gca,s.xlabel{ind})
   ylabel(gca,['p(y|x(',s.xlabel{ind},'))'])
end
end

function [lh,px,px0,X11,X22]=lh2(s,y,x1,x2,ind);
%LH2 computes the two-dimensional log likelihood function over a grid
%   [lh,px,px0,X1,X2]=lh2(s,y,x1,x2,ind);
%
%   ind  {[1 2]}  A vector with two integers indicating which states to vary
%   x1,x2         Vectors with grid points to evaluate lh over
%   y             Measurement(s) as sig object
%   lh            Likelihood function
%   px0           Prior evaluated in x1,x2
%   px            Posterior=px0*lh
%
%   Symbolically, px = log prod_{t=y.t} s.pe(y(t);s.h(t,x,u(t),s.th|x=(x1,x2)))
%
%   Example:
%   s=exsensor('toa',2,1);
%   y=simulate(s,1);
%   plot(s);
%   hold on, lh2(s,y); hold off
%   Larger example:
%   M=5; N=1;
%   s=exsensor('doa',M,1);
%   s.x0=[0.5;0.5];  s.pe=0.1*eye(M);
%   y=simulate(s);
%   subplot(2,2,1),  plot(s),
%   hold on, lh2(s,y);   hold off,   axis([0 1 0 1])
%   subplot(2,2,3),   lh1(s,y,[],1); set(gca,'Xlim',[0 1])
%   subplot(2,2,2),  [p,dum,dum,x]=lh1(s,y,[],2);
%   plot(p,x),set(gca,'Ylim',[0 1])
%   subplot(2,2,4), crlb(s); axis([0 1 0 1])


if nargin<5
   ind=[1 2];
end
if nargin<4
   x2=s.x0(ind(2))+linspace(-1,1,30);
end
if nargin<3
   x1=s.x0(ind(1))+linspace(-1,1,30);
end
if length(ind)~=2
   error('SIGMOD.LH2: ind must be a vector of length 2')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SIGMOD.LH2: ind must be a vector with positive integers')
end
if any(ind)>s.nn(1)
   error('SIGMOD.LH1: state index ind cannot exceed nx in s.nn(1)')
end
if isempty(s.pe)
   error('SIGMOD.LH2: pe cannot be empty for lh2')
end
[X1,X2]=meshgrid(x1,x2);
X11=X1(:);
X22=X2(:);
N1=length(x1);
N2=length(x2);
N=N1*N2;

X=repmat(s.x0',N,1);
X(:,ind)=[X1(:) X2(:)];
if ~isempty(s.px0)
   px0=pdf(s.px0,X);
else
   px0=zeros(N,1);
end
px=px0;
for k=1:length(y.t);
   yhat=feval(s.h,y.t(k),X,y.u(k,:),s.th);
   if size(yhat,2)~=N
      yhat=zeros(s.nn(3),N);
      for n=1:N
         yhat(:,n)=feval(s.h,y.t(k),X(n,:)',y.u(k,:),s.th);
      end
   end
   epsi=repmat(y.y(k,:).',1,N)-yhat;
   lh =pdf(s.pe,epsi.');
   px=px+lh;
%yhat,lh,pause
end
px0=reshape(px0,N1,N2);
px=reshape(px,N1,N2);
lh=reshape(lh,N1,N2);
if nargout==0
%   surf(lh)
   contour2(x1,x2,log(lh))
%,'linewidth',2)
   xlabel(gca,s.xlabel{ind(1)})
   ylabel(gca',s.xlabel{ind(2)})
end
end

% Plot tools
% -----------

function plot(varargin)
%PLOT illustrates a sensor network represented by sigmod s
%   plot(s1,s2,...,'Property','Value')
%
%   The sensor locations are by default assumed to be in the field s.th,
%   whose size is twice the number of sensors.
%   The target locations are by default assumed to be in the field s.x,
%   whose size is twice the number of sensors.
%   Both these assumptions can be changed by an index property.
%
%   Property  Value Description
%   -------------------------------------------------------------
%   xind      {:}   Target position indeces. Use [1 2] for one target with
%                   many states
%   thind     {:}   Sensor position indeces. Use [1 2] for one sensor when
%                   there are more parameters in th
%
%
%   Example:
%   s=exsensor(1,5);
%   plot(s)


N=0;  K=0; optvar='';
k=0;
while k<length(varargin)
    k=k+1;
    if isstr(varargin{k})  % Property Value Pair
        K=K+1;
        optvar{K}=varargin{k};
        K=K+1;
        k=k+1;
        optvar{K}=varargin{k};
    else   % vector or struct
        N=N+1;
        Z{N}=varargin{k};
    end
end
opt=struct('axis',0,'scatter','off','conf',0,'Xlim',[],'Ylim',[],'col','bgrmyk','linewidth',[],'fontsize',[],'view','staircase','legend',[],'ind',':','xind',[],'thind',[]);

opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

Ncol=length(opt.col);
if N>Ncol  % repeat the color pattern
   opt.col=char(kron(ones(1,ceil(N/Ncol)),opt.col));
end

for n=1:N;
   s=Z{n};
   if ~( isa(s,'sigmod') | isa(s,'nl'))
      error('SIGMOD.PLOT: input argument must be a NL/SIGMOD object')
   end
   if isempty(opt.thind) | opt.thind==':'
      p=s.th(:);
      Pp=s.P;
   else
      p=s.th(opt.thind);
      Pp=s.P(opt.thind,opt.thind);
   end
   ns=length(p)/2;
   p=reshape(p,2,ns)';
   if isempty(opt.xind) | opt.xind==':'
      x=s.x0;
   else
      x=s.x0(opt.xind);
   end
   nt=length(x)/2;
   x=reshape(x,2,nt)';
   for k=1:nt
      plot(x(k,1),x(k,2),['*',opt.col(n)],'linewidth',opt.linewidth)
      hold on
      text(x(k,1),x(k,2),['T',num2str(k)],'color',opt.col(n),'fontsize',opt.fontsize)
      if ~isempty(s.px0)
           P=cov(s.px0);
	   plot2(ndist(x(k,:)',1e14*P(2*k-1:2*k,2*k-1:2*k)),'col',opt.col(n),'levels',1,'legend','off')
      end
   end
   for k=1:ns
	   plot(p(k,1),p(k,2),['*',opt.col(n)],'color',opt.col(n))
           if ~isempty(Pp)
	      plot2(ndist(p(k,:)',Pp(2*k-1:2*k,2*k-1:2*k)),'col',opt.col(n),'levels',1,'legend','off')
           end
	   text(p(k,1),p(k,2),['S',num2str(k)],'color',opt.col(n))
   end
end
axis([ 0.8*min([p(:,1);x(:,1)]) 1.25*max([p(:,1);x(:,1)]) 0.8*min([p(:,2);x(:,2)]) 1.25*max([p(:,2);x(:,2)]) ])
if ~oldhold
   hold off
end
end

end %methods
end %sigmod
