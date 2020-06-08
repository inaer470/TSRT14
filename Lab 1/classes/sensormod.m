classdef sensormod < nl

%SENSORMOD is the model object for sensor models
%   sensormod is a special case of the NL object with f='x'
%   m=sensormod(h,nn)
%   The signal model is defined as
%        y(t) = h(t,x(t),u(t);th) + e(t),  e(t)~pe
%     E[x(0)] = x0, x(0)~px0
%   h is an inline object or a string with arguments t,x,u,th,
%   and nn=[nx nu ny nth] gives the dimensions of x,u,y and the parameters th.
%   Convention of dimensions:
%     x (nx,1) vector or (nx,N) matrix
%     u (nu,1) vector or (nu,N) matrix
%     y (ny,1) vector or (ny,N) matrix
%     th vector of length nth
%
%   Distributions pe, px0 are specified as pdfclass objects, or
%   as a covariance matrix of suitable dimension, in which
%   case a ndist object is created.
%
%   help sensormod.sensormod    gives help for the constructor
%   methods sensormod           lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $


methods

function m=sensormod(h,nn)
%SENSORMOD constructor signal models
%   m=sensormod(h,nn)
%   The signal model is defined as
%        y(t) = h(t,x(t),u(t);th) + e(t),  e(t)~pe
%     E[x(0)] = x0, x(0)~px0
%   h is an inline object or a string with arguments t,x,u,th,
%   and nn=[nx nu ny nth] gives the dimensions of x,u,y and the parameters th.
%   Convention of dimensions:
%     x (nx,1) vector or (nx,N) matrix
%     u (nu,1) vector or (nu,N) matrix
%     y (ny,1) vector or (ny,N) matrix
%     th vector of length nth
%
%   Distributions pe, px0 are specified as pdfclass objects, or
%   as a covariance matrix of suitable dimension, in which
%   case a ndist object is created.
%   Other fields: name, xlabel thlabel ulabel ylabel
%
%   Examples:
%   s=exsensor(1,5);
%   s=sensormod('sin(2*pi*th(1)*t)',[0 0 1 1]);
%   s.pe=0.01;  s.th=0.1;
%   s.thlabel='Frequency'; s.name='Sine wave in noise';
%   y=simulate(s,1:10);
%   shat=estimate(s,y)

%#
%#

if isa(h,'sensormod')
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
   error('sensormod: two input arguments required')
else
   mnl=nl('x',h,nn,NaN);
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
%SIMULATE simulates a sensor model
%
%   1. y=simulate(s,t), gives z=s(t,x) at times t and state s.x0.
%   2. y=simulate(s,x), gives z=s(t,x) at times x.t and state x.x.
%
%   s sensormod object
%   t time instant or vector (default 0), gives y=s.h(t,s.x0)
%   t = x sig object, gives y=s.h(x.t,x.x)
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
      error('SENSORMOD.simulate: time argument required, s(0) not defined')
   end
end
opt=struct('MC',30);
opt=optset(opt,varargin);
if isa(t,'sig')
  tk=t.t;
  N=length(tk);
  u=t.u;
  x=t.x;
else
  tk=t;
  N=length(tk);
  u=zeros(N,0);
  x=repmat(s.x0(:)',N,1);
end

y=[];
ny=s.nn(3);
y=zeros(ny,N);
try
  yk=feval(s.h,tk,x,u,s.th);
  yk=reshape(ny,N/ny);
catch
  for k=1:N
    yk=feval(s.h,tk(k),x(k,:)',u(k,:)',s.th);
    y(:,k)=yk(:);
  end
end
if isa(s.pe,'pdfclass')
  y=y+rand(s.pe,N)';
end

yMC=[];
xMC=[];
if ~isempty(s.P) & any(any(s.P~=0))
       nx=s.nn(1); nu=s.nn(2); ny=s.nn(3); nth=s.nn(4); N=length(t);
       yMC=zeros(opt.MC,N,ny);
       xMC=zeros(opt.MC,N,nx);
       thMC=rand(ndist(zeros(nth,1),s.P(1:nth,1:nth)),opt.MC)';
       mtmp=s;
       for i=1:opt.MC;
           mtmp.th=s.th+thMC(:,i);
           ztmp=simulate(mtmp,t,varargin{:},'MC',0);
           yMC(i,:,:)=ztmp.y;
           xMC(i,:,:)=ztmp.x;
       end
end
z=sig(y.',tk,u,x,yMC,xMC);
z.fs=s.fs;
z.MC=opt.MC;
z.xlabel=s.xlabel;
z.ylabel=s.ylabel;
z.ulabel=s.ulabel;
end

%==================
%----Detectors-----
%==================

function [b,level,h,T]=detect(s,y,pfa)
%DETECT evaluates a hypothesis test H0: y~e vs H1: y~s+e
%   [b,level,h,T]=detect(s,y,pfa)
%
%      H0: y=e,
%      H1: y=s.h(x0)+e,
%   where e~N(0,R), R=cov(s.pe)
%   The test statistic is distributed as
%      H0: T(y) = y'*inv(R)*y ~ chi2dist(nx)
%      H1: T(y) = y'*inv(R)*y ~ ncchi2dist(nx,lambda)
%   where lambda = s.h(x0)'*inv(R)*s.h(x0)
%   Detection formula: Pd=1-erf(erfinv(1-Pfa)-lambda)
%
%   y     data sample as SIG model or vector
%   pfa   false alarm rate (default 0.01)
%   b     binary decision
%   level level of the test
%   h     threshold corresponding to pfa
%   T     Test statistic T(y)=y'*inv(R)*y
%
%   Example:
%   s=exsensor('toa',5,1);
%   y=simulate(s)
%   [b,l,h,T]=detect(s,y,0.01);
%
%   See also: lrt, roc, pdplot1, pdplot2

if nargin<3, pfa=0.01; end
if nargin<2, error('SENSORMOD.DETECT: a third data argument y required'), end

Z=chi2dist(s.nn(1));
h=erfinv(Z,1-pfa);
T=y.y*inv(cov(s.pe))*y.y';
b=(T>h);
level=erf(Z,T);

end

function [pd,pfa]=roc(s,x0,h)
%ROC plots the receiver operating characteristics curve Pd as a function of Pfa
%   [pd,pfa]=roc(s,x0,h)
%   Test the hypotheses
%   The hypothesis test is
%      H0: y=e,
%      H1: y=s.h(x0)+e,
%   where e~N(0,R), R=cov(s.pe)
%   The test statistic is distributed as
%      H0: T(y) = y'*inv(R)*y ~ chi2dist(nx)
%      H1: T(y) = y'*inv(R)*y ~ ncchi2dist(nx,lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%   Detection formula: Pd(Pfa)=1-erf(erfinv(1-Pfa)-lambda)
%
%   x0   is a nx times nroc matrix, where nroc is the number of ROC curves
%        default (if x0 omitted or empty) x0=s.x0
%   h    is a vector of threshold. Default it is gridded as h(pfa)
%
%   With two input arguments, h is gridded automatically
%   Without output arguments, a plot is generated
%
%   Example:
%    s=exsensor('toa',5,1); % Default network
%    s.pe=1*eye(5);         % Increase noise level
%    roc(s,[1 1;0.5 0.5]'); % Two ROC curves
%
%   See also: detect, pdplot2

if nargin<2 | isempty(x0)
   x0=s.x0;
end
[nx,nroc]=size(x0);
if nx~=s.nn(1)
   error('SENSORMOD.ROC: x0 must be a nx times nroc matrix')
end
Z=chi2dist(s.nn(1));
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

for j=1:nroc
   k=1; %time
%   H=numgrad(s.h,2,0,x0(:,j),[],s.th)
%   T=x0(:,j)'*H'*inv(cov(s.pe))*H*x0(:,j);
   hx0=feval(s.h,0,x0(:,j),[],s.th);
   T=hx0'*inv(cov(s.pe))*hx0;
   for i=1:N
      Z1=ncchi2dist(s.nn(1),T);
      try
         pd(i,j)=1-erf(Z1,h(i));
      catch
         pd(i,j)=1;
      end
   end
end
if nargout==0
   lw=2;
   fs=14;
   ind=1:round(N/5):N;
   plot(pfa,pd,'-',pfa(ind),pd(ind,:),'.','linewidth',lw)
   hold on
   for j=1:nroc
      for i=1:length(ind);
         text(pfa(ind(i)),pd(ind(i),j),num2str(h(ind(i)),'%.3g'),'fontsize',fs)
      end
   end
   hold off
   ylabel(gca,'Pd')
   xlabel(gca,'Pfa')
end
end

function pd=pdplot1(s,x1,pfa,ind)
%PDPLOT1 plots detection probability Pd as a function on the grid x1 for given Pfa
%   pd=pdplot1(s,x1,pfa,ind)
%
%   The hypothesis test is
%      H0: y=e,
%      H1: y=s.h(x0)+e,
%   where e~N(0,R), R=cov(s.pe)
%   The test statistic is distributed as
%      H0: T(y) = y'*inv(R)*y ~ chi2dist(nx)
%      H1: T(y) = y'*inv(R)*y ~ ncchi2dist(nx,lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%
%   For both models, x0(ind) is gridded according to x1
%   The output gives pd from the formula:
%      Pd=1-erf(erfinv(1-Pfa)-lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%   Without output argument, a plot of pd(x1) is generated.
%
%   x1    grid vector
%   pfa   False alarm rate
%   ind   Index of s.x0 to vary, that is, s.x0(ind)=x1(i)
%
%   Example:
%    s=exsensor('toa',5,1);       % Default network
%    s.pe=1*eye(5);               % Increase noise level
%    pdplot1(s,0:0.05:1);      % Plot
%
%   See also: pdplot2, detect, roc

if nargin<4
   ind=[1];
end
if nargin<3
   pfa=0.01;
end
if pfa<=0 | pfa >=1
   error('sensormod.pdplot1: pfa must be in (0,1)')
end
if nargin<2
   error('sensormod.pdplot1: Two input arguments needed at least, pdplot1(s,x1)')
end

x0=s.x0;

for i=1:length(x1)
      k=1; %time
      Z0=chi2dist(s.nn(1));
      x0(ind(1))=x1(i);
         %H0=numgrad(s0.h,2,k,x00,[],s0.th);
         %T0=x00'*H0'*inv(cov(s0.pe))*H0*x00;
      hx=feval(s.h,0,x0,[],s.th);
      T=hx'*inv(cov(s.pe))*hx;
      Z1=ncchi2dist(s.nn(1),T);

      h=erfinv(Z0,1-pfa);
      try
         pd(i)=1-erf(Z1,h);
      catch
         pd(i)=1;
      end
      if isnan(pd(i)), pd(i)=0; end
end
if nargout==0
   plot(x1,pd,'-',s.x0(ind(1)),0,'d','markersize',16)
   if ~isempty(s.name)
      title(gca,[s.name,', Pfa = ',num2str(pfa)])
   else
      title(gca,['Pfa = ',num2str(pfa)])
   end
   if ~isempty(s.xlabel)
      ylabel(gca,'P_D')
      xlabel(gca,s.xlabel{ind(1)})
   end
end
end

function pd=pdplot2(s,x1,x2,pfa,ind)
%PDPLOT2 plots detection probability Pd as a function on the grid x1,x2 for given Pfa
%   pd=pdplot2(s,x1,x2,pfa,ind)
%
%   The hypothesis test is
%      H0: y=e,
%      H1: y=s.h(x0)+e,
%   where e~N(0,R), R=cov(s.pe)
%   The test statistic is distributed as
%      H0: T(y) = y'*inv(R)*y ~ chi2dist(nx)
%      H1: T(y) = y'*inv(R)*y ~ ncchi2dist(nx,lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%
%   For both models, x0(ind) is gridded according to x1, x2
%   The output gives pd from the formula:
%      Pd=1-erf(erfinv(1-Pfa)-lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%   Without output argument, a contour plot of pd(x1,x2) is generated.
%
%   x1    grid vector
%   pfa   False alarm rate
%   ind   Index of s.x0 to vary, that is, s.x0(ind)=x1(i)
%
%   Example:
%    s=exsensor('toa',5,1);       % Default network
%    s.pe=1*eye(5);               % Increase noise level
%    pdplot2(s,0:0.05:1,0:0.05:1);      % Plot
%
%   See also: pdplot2, detect, roc

if nargin<5
   ind=[1 2];
end
if nargin<4
   pfa=0.01;
end
if pfa<=0 | pfa >=1
   error('sensormod.pdplot2: pfa must be in (0,1)')
end
if nargin<3
   error('sensormod.pdplot2: Three input arguments needed at least, pdplot2(s,x1,x2)')
end

x0=s.x0;

for i=1:length(x1)
   for j=1:length(x2)
      k=1; %time
      Z0=chi2dist(s.nn(1));
      x0(ind)=[x1(i) x2(j)];
      hx=feval(s.h,0,x0,[],s.th);
      T=hx'*inv(cov(s.pe))*hx;
      Z1=ncchi2dist(s.nn(1),T);
      h=erfinv(Z0,1-pfa);
      try
         pd(i,j)=1-erf(Z1,h);
      catch
         pd(i,j)=1;
      end
      if isnan(pd(i,j)), pd(i,j)=0; end
   end
end
if nargout==0
   lev=1-logspace(-1,0,10);
   [cs,h]=contour(x1,x2,pd); %,lev,'linewidth',lw)
   clabel(cs,h)
   hold on
   plot(s.x0(ind(1)),s.x0(ind(2)),'d','markersize',16)
   hold off
   if ~isempty(s.name)
      title(gca,[s.name,', Pfa = ',num2str(pfa)])
   else
      title(gca,['Pfa = ',num2str(pfa)])
   end
   if ~isempty(s.xlabel)
      ylabel(gca,s.xlabel{ind(2)})
      xlabel(gca,s.xlabel{ind(1)})
   end
end
end


function pd=pdplot3(s1,s0,x1,x2,pfa,ind)
%Obsolute function
% PDPLOT2 plots detection probability Pd as a function on the grid x1, x2 for given Pfa
%   pd=pdplot2(s1,s0,x1,x2,pfa,ind)
%
%   The hypothesis test is
%      H0: y=s0.h(s0.x0)+e0, e0~N(0,R0), R0=cov(s0.pe)
%      H1: y=s1.h(s1.x0)+e1, e1~N(0,R1), R1=cov(s1.pe)
%   Default s0.h=0 and R0=R1 if s0=[]
%   For both models, x0(ind) is gridded according to x1,x2
%   The output gives pd from the formula:
%      Pd=1-erf(erfinv(1-Pfa)-lambda)
%   where lambda=x0'*H'*inv(R)*H*x0
%   Without output argument, a contour plot pd(x1,x2) is generated.
%
%   ind   Indexes of s.x0 to vary, that is, s.x0(ind)=[x1(i),x2(j)];
%   pfa   False alarm rate
%   ind   Indexes s.x0 to the variable states
%
%   Example:
%    s=exsensor('toa',5,1);           % Default network
%    s.pe=1*eye(5);                   % Increase noise level
%    pdplot2(s,[],0:0.05:1,0:0.05:1); % Contour plot
%
%   See also: pdplot1, detect, roc

if nargin<6
   ind=[1 2];
end
if nargin<5
   pfa=0.01;
end

x01=s1.x0;
if ~isempty(s0)
    x00=s0.x0;
end
for i=1:length(x1)
   for j=1:length(x2)
      k=1; %time
      if isempty(s0)
         Z0=chi2dist(s1.nn(1));
      else
         x00(ind(1))=x1(i);
         x00(ind(2))=x2(j);
         H0=numgrad(s0.h,2,k,x00,[],s0.th);
         T0=x00'*H0'*inv(cov(s0.pe))*H0*x00;
         Z0=ncchi2dist(s0.nn(1),T0);
      end
      h=erfinv(Z0,1-pfa);
      x01(ind(1))=x1(i);
      x01(ind(2))=x2(j);
      H1=numgrad(s1.h,2,k,x01,[],s1.th);
      T1=x01'*H1'*inv(cov(s1.pe))*H1*x01;
      Z1=ncchi2dist(s1.nn(1),T1);
      try
         pd(i,j)=1-erf(Z1,h);
      catch
         pd(i,j)=1;
      end
  %    Z0,Z1,T0,T1,H0,H1,x00,x01,h,pd(i,j),pause
      if isnan(pd(i,j)), pd(i,j)=0; end
   end
end
if nargout==0
   lev=1-logspace(-1,0,10);
   [cs,h]=contour(x1,x2,pd); %,lev,'linewidth',lw)
   clabel(cs,h)
   if ~isempty(s1.name)
      title(gca,[s1.name,', Pfa = ',num2str(pfa)])
   else
      title(gca,['Pfa = ',num2str(pfa)])
   end
   if ~isempty(s1.xlabel)
      ylabel(gca,s1.xlabel{ind(2)})
      xlabel(gca,s1.xlabel{ind(1)})
   end
end
end


%==================
%----Estimators----
%==================


function [xhat,shat]=ls(s,y);
%LS computes the least squares estimate
%   [xhat,shat]=ls(s,y);
%
%   The function computes
%   argmin V(x)= argmin sum( (y-H*x)^T*(y-H*x))
%   For nonlinear functions s.h, numgrad is used to
%   linearize around s.x0, which then becomes a crucial parameter.
%
%   Example:
%   s=exsensor('toa',3,1);
%   y=simulate(s,1);
%   xhat=ls(s,y);
%   plot(s)
%   hold on
%   xplot2(xhat,'conf',90)
%   hold off

if isempty(s.pe)
   error('sensormod.ls: pe must be defined')
end
for k=1:length(y.t)
   H=numgrad(s.h,2,y.t(k),s.x0,y.u(k,:)',s.th);
   yhat=bsxfun(@plus, feval(s.h,y.t(k),s.x0,y.u(k,:)',s.th), mean(s.pe));
   PP=pinv(H'*H);
   ff=H'*(y.y(k,:)'-yhat+H*s.x0);
   x(k,:) = PP*ff;
   yy(k,:) = yhat'+x(k,:)*H'-s.x0'*H';
   yy(k,:)=yhat;  %+x(k,:)*H';
   PP=pinv(H'*H);
   Px(k,:,:)=PP*(H'*cov(s.pe)*H)*PP;
   Py(k,:,:)=H*squeeze(Px(k,:,:))*H';
end
xhat=sig(yy,y.t,y.u,x,Py,Px);
if length(y.t)==1
   shat=s;
   shat.x0=x';
   shat.px0=squeeze(Px(1,:,:));
end
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
%   s=exsensor('toa',3,1)
%   y=simulate(s,1)
%   [xhat,shat]=wls(s,y); shat  % Note x0 distribution
%   plot(s)
%   hold on
%   plot(shat)                  % overlaid sensor network
%   y=simulate(s,3);
%   [xhat,shat]=wls(s,y);       % Three estimates
%   xplot2(xhat,'conf',90)
%   hold off

if isempty(s.pe)
   error('sensormod.wls: pe must be defined')
end
R=cov(s.pe);
Rinv=pinv(R);
for k=1:length(y.t)
   H=numgrad(s.h,2,y.t(k),s.x0,y.u(k,:)',s.th);
   yhat=bsxfun(@plus, feval(s.h,y.t(k),s.x0,y.u(k,:)',s.th), mean(s.pe));
   PP=pinv(H'*Rinv*H);
   ff=H'*Rinv*(y.y(k,:)'-yhat+H*s.x0);
   x(k,:) = PP*ff;
   yy(k,:) = yhat'+x(k,:)*H'-s.x0'*H';
   PP=0.5*(PP+PP');
   Px(k,:,:)=PP;
   Py(k,:,:)=H*PP*H';
end
xhat=sig(yy,y.t,y.u,x,Py,Px);
if length(y.t)==1
   shat=s;
   shat.x0=x';
   shat.px0=squeeze(Px(end,:,:));
else
    shat=[];
end
end

function [xhat,shat,res]=ml(s,y);
%ML computes the ML/NLS parameter estimate in the sensormod object s from y
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
try
   yhat.x=shat.x0';
end
xhat=crlb(nl(shat),yhat);
%xhat=sig(yhat.y,y.t,y.u,yhat.x);
end

function x=crlb(s,y)
%CRLB computes the Cramer-Rao lower bound for the state at y.x
%   x=crlb(s,y)
%
%   s   is a sensormod object
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
if nargin<2 | isempty(y)
    x0=s.x0';
    u=[];
else
  if length(y)>1
     error('CRLB can only be computed for one sample at the time')
  end
  x0=y.x;
  u=y.u;
end
y=sig(x0,s.fs,u,x0);
x=crlb(nl(s),y);
if nargout==0
   xplot2(x,'conf',90)
end
end

function [cx,X1,X2]=crlb2(s,y,x1,x2,ind,type);
%CRLB2 computes a scalar CRLB measure over a 2D state space grid
%   [cx,X1,X2]=crlb2(s,y,x1,x2,ind,type);
%
%   s   is a sensormod object
%   y   is a sig object of length 1 specifying x0 as y.x.
%       If y is empty, x0 is taken from s.x0
%   x1,x2 contain the grid for the states x(ind) (values in x0 are
%   omitted for these two dimensions)
%   ind  {[1 2]}  A vector with two integers indicating which
%                 states x1,x2 refer to in x
%   type  'trace', {'rmse'}, 'det', 'max'
%                 operator for transforming P(x) to scalar c(x)
%                 rmse=sqrt(trace)
%
%   Symbolically, cx=type(crlb(s,x))
%
%   Example:
%   s=exsensor('toa',2,1);
%   s.pe=1*eye(2);
%   plot(s);
%   hold on, crlb2(s); hold off

if nargin<6
   type='rmse';
end
if nargin<5
   ind=[1 2];
end
if nargin<4
   x2=linspace(0,1,22);
end
if nargin<3
   x1=linspace(0,1,20);
end
if nargin<2
   y=[];
end
if length(ind)~=2
   error('SENSORMOD.CRLB2: ind must be a vector of length 2')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SENSORMOD.CRLB2: ind must be a vector with positive integers')
end
if  any(ind>s.nn(1))
   error('SENSORMOD.CRLB2: The state indeces in ind cannot exceed nx=s.nn(1)')
end
if isempty(s.pe)
   error('SENSORMOD.CRLB2: pe cannot be empty for crlb2')
end
if ~isstr(type)
   error('SENSORMOD.CRLB2: type must be a string')
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
     error('SENSORMOD.CRLB2 can only be computed for one sample at the time')
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
       error(['SENSORMOD.CRLB2: unknown option type = ',type])
   end
end
cx=reshape(cx,N2,N1);
if nargout==0
    cxsort=sort(cx(:));
    lev=cxsort(10:N/10:N);
   [cs,h]=contour(x1,x2,cx,lev);
   clabel(cs,h)
   xlabel(gca,s.xlabel{ind(1)})
   ylabel(gca,s.xlabel{ind(2)})
end
end

function [lh,x1,px,px0]=lh1(s,y,x1,ind);
%LH1 computes the one-dimensional likelihood function over a state space grid
%   [lh,x1,px,px0]=lh1(s,y,x,ind);
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
   error('SENSORMOD.LH1: ind must be an integer')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SENSORMOD.LH1: ind must be a positive integer')
end
if ind>s.nn(1)
   error('SENSORMOD.LH1: state index ind cannot exceed nx in s.nn(1)')
end
if isempty(s.pe)
   error('SENSORMOD.LH1: pe cannot be empty for lh1')
end
if s.nn(1)~=y.nn(1)
   error(['SENSORMOD.LH1: model and data do not have the same ' ...
          'measurement dimension ny'])
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
   ylabel(gca,['p(y|',s.xlabel{ind},')'])
end
end

function [lh,x1,x2,px,px0,X1,X2]=lh2(s,y,x1,x2,ind);
%LH2 computes the two-dimensional log likelihood function over a grid
%   [lh,x1,x2,px,px0,X1,X2]=lh2(s,y,x1,x2,ind);
%
%   Symbolically, the notation is defined as
%     x=x0; x0(ind)=[x1 x2]; lh(x)=p(x|y)
%
%   ind  {[1 2]}  A vector with two integers indicating which states to vary
%   x1,x2         Vectors with grid points to evaluate lh over
%   y             Measurement(s) as sig object
%   lh            Likelihood function evaluated in (x1,x2)
%   px0           Prior evaluated in (x1,x2)
%   px            Posterior=px0*lh
%   X1,X2         Meshgrid
%
%   The likelihood is computed as
%      px = log prod_{t=y.t} s.pe(y(t);s.h(t,x,u(t),s.th|x=(x1,x2)))
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
%   s.th=[0.1 0.5 0.6 0.9 0.6 0.1 0.2 0.8 0.2 0.2];
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
   x2=s.x0(ind(2))+linspace(-0.5,0.5,30);
end
if nargin<3
   x1=s.x0(ind(1))+linspace(-0.5,0.5,30);
end
if s.nn(3)~=y.nn(3)
   error(['SENSORMOD.LH2: model and data do not have the same ' ...
          'measurement dimension ny'])
end
if length(ind)~=2
   error('SENSORMOD.LH2: ind must be a vector of length 2')
end
if ~isequal(ind,round(ind)) | any(ind<=0)
   error('SENSORMOD.LH2: ind must be a vector with positive integers')
end
if any(ind)>s.nn(1)
   error('SENSORMOD.LH2: state index ind cannot exceed nx in s.nn(1)')
end
if isempty(s.pe)
   error('SENSORMOD.LH2: pe cannot be empty for lh2')
end
[X1, X2]=meshgrid(x1, x2);
sz = size(X1);
N = numel(X1);

X =  repmat(s.x0', [N, 1]);
X(:,ind)=[X1(:) X2(:)];
if ~isempty(s.px0)
   px0=pdf(s.px0,X);
else
   px0=zeros(N,1);
end
px=px0;
for k=1:length(y.t);
    yhat=[];
    try
        yhat=feval(s.hh,y.t(k),X,y.u(k,:),s.th);
    end
   if size(yhat,2)~=N
      yhat=zeros(s.nn(3),N);
      for n=1:N
         yhat(:,n)=feval(s.h,y.t(k),X(n,:)',y.u(k,:),s.th);
      end
   end
   epsi=repmat(y.y(k,:).',1,N)-yhat;
   %   lh =-log(pdf(s.pe,epsi.')); % Numerically unstable
   for n=1:N;
        lh(n)=epsi(:,n)'*inv(cov(s.pe))*epsi(:,n); % assumes Gaussian pe
   end
   px=px+lh(:);
end
lhind=find(lh==Inf);
lh(lhind)=10*max(lh(:));
px0=reshape(px0, sz);
px=reshape(px, sz);
lh=reshape(lh, sz);

if nargout==0
    %v=linspace(min(lh(:)), median(lh(:)),20);
    vs=sort(lh(:));
    v=vs(round(N/30):round(N/30):round(N/5));
   contour2(x1,x2,lh,v)
   xlabel(gca,s.xlabel{ind(1)})
   ylabel(gca',s.xlabel{ind(2)})
end
end

% Plot tools
% -----------

function plot(varargin)
%PLOT illustrates a sensor network represented by sensormod s
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
%   conf      {90}  Confidence interval for targets AND sensors
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
opt=struct('axis',0,'scatter','off','conf',0,'Xlim',[],'Ylim',[],'col','bgrmyk','linewidth',[],'fontsize',[],'markersize',[],'view','staircase','legend',[],'ind',':','xind',':','thind',[]);

opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

Ncol=length(opt.col);
if N>Ncol  % repeat the color pattern
   opt.col=char(kron(ones(1,ceil(N/Ncol)),opt.col));
end
x1=[]; x2=[];
for n=1:N;
   s=Z{n};
   if ~( isa(s,'sensormod') | isa(s,'nl'))
      error('SENSORMOD.PLOT: input argument must be a NL/SENSORMOD object')
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
   nt=floor(length(x)/2);
   x=reshape(x(1:2*nt),2,nt)';
   hold on;
   box on;
   for k=1:nt
      plot(x(k,1),x(k,2),['*',opt.col(n)],'linewidth',opt.linewidth)
      text(x(k,1),x(k,2),['T',num2str(k)],'color',opt.col(n),'fontsize',opt.fontsize)
      if ~isempty(s.px0)  & opt.conf>0
           P=cov(s.px0);
	   plot2(ndist(x(k,:)',P(2*k-1:2*k,2*k-1:2*k)),'col',opt.col(n),'conf',opt.conf,'legend','off','linewidth',opt.linewidth);
      end
   end
   for k=1:ns
	   plot(p(k,1),p(k,2),['*',opt.col(n)],'color',opt.col(n),'markersize',opt.markersize)
           if ~isempty(Pp) & opt.conf>0
	      plot2(ndist(p(k,:)',Pp(2*k-1:2*k,2*k-1:2*k)),'col',opt.col(n),'conf',opt.conf,'legend','off','linewidth',opt.linewidth);
           end
	   text(p(k,1),p(k,2),['S',num2str(k)],'color',opt.col(n),'fontsize',opt.fontsize)
   end
   x1=[p(:,1);x(:,1)];
   x2=[p(:,2);x(:,2)];
end
set(gca,'linewidth',opt.linewidth,'fontsize',opt.fontsize)
if all(x1==0)
    x1lim=[-1 1];
else
    x1lim=[(0.8+0.225*(1-sign(min(x1))))*min(x1) (0.8+0.225*(1+sign(max(x1))))*max(x1)];
end
if all(x2==0)
    x2lim=[-1 1];
else
    x2lim=[(0.8+0.225*(1-sign(min(x2))))*min(x2) (0.8+0.225*(1+sign(max(x2))))*max(x2)];
end
axis([x1lim x2lim])

if ~oldhold
   hold off
end
end

end %methods
end %sensormod
