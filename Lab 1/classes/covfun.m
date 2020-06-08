classdef covfun
%COVFUN is the covariance object
%   Fields:
%     R    is (taumax+1,ny,ny)
%     tau  is a vector
%     RMC  is (MC,taumax+1,ny,ny)
%   COVFUN constructor
%     R=covfun              % empty object
%     R=covfun(R,tau)       % R(tau,:,:)
%     R=covfun(R,tau,RMC)   % with MC simulations
%     R=covfun(s)           % Conversions from SIG and LTI (ARX, SS) objects
%   Optional public fields: fs MC name xlabel ulabel ylabel desc
%   Change these with R.fs=fs etc.
%   Examples:

% Copyright Fredrik Gustafsson, Sigmoid AB
% $ Release: $

properties (SetAccess = public)
   R, RMC, tau
   fs, MC, name='empty'; tlabel, ylabel, ulabel, xlabel, desc
   marker, markerlabel
end

methods

function R=covfun(varargin)
%COVFUN constructor
%   cov=covfun              % empty object
%   cov=covfun(R,tau)       % R(tau,:,:)
%   cov=covfun(R,tau,RMC)   % with MC simulations
%   cov=covfun(s)           % Conversions from SIG and LTI (ARX, SS) objects
%   Optional public fields: fs MC name xlabel ulabel ylabel desc
%   Change these with R.fs=fs etc.
%
%   Example:
%     N=1000;
%     y=filter(1,[1 1.6 0.64],randn(N,1));
%     for tau=0:10;
%         R(1+tau)=sum(y(1:N-tau)'*y(1+tau:N))/(N-tau);
%     end
%     c1=covfun(R,0:10);
%     c1.name='AR(2) example';
%     c1.tlabel='Time [samples]';
%     c1.ylabel='Volt';
%     ysig=sig(y);
%     c2=covfun(ysig);
%     subplot(2,1,1), plot(c1)
%     subplot(2,1,2), plot(c2)

%#
%#

R=[]; tau=NaN; RMC=[];
fs=NaN; MC=0;
name='';  xlabel=[]; ulabel=[]; ylabel=[]; desc=[];
marker=[]; makerlabel='';

if nargin<1
    % Do nothing
elseif isa(varargin{1},'covfun')  % Copy contents
    c=varargin{1};
    R=c.R; tau=c.tau; fs=c.fs; name=c.name; RMC=c.RMC;
elseif isa(varargin{1},'lss')
    c=ss2covfun(varargin{1},varargin{2:end});
    R=c.R; tau=c.tau; fs=c.fs; name=c.name; RMC=c.RMC;
elseif isa(varargin{1},'arx')
    c=arx2covfun(varargin{1},varargin{2:end});
    R=c.R; tau=c.tau; fs=c.fs; name=c.name; RMC=c.RMC;
elseif isa(varargin{1},'sig')
    c=sig2covfun(varargin{1},varargin{2:end});
    R=c.R;  tau=c.tau;fs=c.fs; name=c.name; RMC=c.RMC;
elseif isnumeric(varargin{1})
    R=varargin{1};
    if ~isnumeric(R), error('R must be numeric'), end
    if size(R,2)~=size(R,3)
        if size(R,3)==1 & size(R,1)==1  % scalar signal, transpose R
            R=R';
        else
            error('R must have dimension (taumax+1,ny,ny)')
        end
    end
    if nargin>1
        tau=varargin{2};
        if ~isnumeric(tau), error('tau must be numeric'), end
        if length(tau)~=size(R,1),
            error('Length of tau must match the size of the first dimension of R')
        end
    end
    if nargin>2
        RMC=varargin{3};
        if ~isempty(RMC)
          if size(RMC,2)~=size(R,1) | size(RMC,3)~=size(R,2) | size(RMC,4)~=size(R,3)
            error('Last three dimensions of RMC must have the same size as R')
          end
        end
        MC=size(RMC,1);
    end
else
   error('Incorrect syntax to COVFUN')
end
%     R,tau,fs,name,RMC
cov.R=R; cov.tau=tau; cov.fs=fs; cov.name=name; cov.RMC=RMC;
end


function disp(R)
disp(['COVFUN object: ',name])
end

function c=estimate(c,s)
%ESTIMATE estimates a covariance function from data
%   The sig method sig2covfun is used
%   Example:
%     m0=rand(arx(4));
%     y=simulate(m0,100);
%     csighat=estimate(covfun,y);  % Equivalent to  csighat=covfun(y);
%     c0=covfun(m0);
%     plot(c0,csighat)
%   See also: sig2covfun
c=covfun(s);
end


function [ny,ntau]=size(R)
%SIZE returns the sizes of R
%   [ny,ntau]=size(R)
ny=size(R.R,2);
ntau=size(R.R,1);
end

%==============================================
%-----OPERATORS FOR STOCHASTIC OBJECTS----------
%==============================================

function Cr=rand(C,n)
%RAND returns a cell array of COVFUN realizations taken from MC data
%   Cr=rand(C,n);
%   n is the number of realizations (default n=1)
%   Cr is a cell array of length n

if nargin<2; n=1; end
if isempty(C.RMC)
   error('No MC data in input argument')
end
MC=size(C.RMC,1);
for k=1:n
   ind=ceil(MC*rand(1,1));
   Cr{k}=covfun(shiftdim(C.RMC(ind,:,:),1),C.tau);
   Cr{k}=inherit(Cr{k},C);
end

end

function C2=fix(C1)
%FIX removes the Monte Carlo data from the COVFUN object
C2=C1;
C2.RMC=[];
C2.MC=0;
end

function C2=mean(C1)
%MEAN returns the mean signal of the Monte Carlo data
C2=C1;
if isempty(C1.RMC)
   error('No MC data in input argument')
end
C2.R=shiftdim(super.mean(C2.RMC,1),1);
C2.RMC=[];
C2.MC=0;
end

function C2=E(C1)
%E returns the mean signal of the Monte Carlo data
C2=mean(C1);
end

function C2=std(C1)
%STD returns the standard deviation of the Monte Carlo data
C2=C1;
if isempty(C1.RMC)
   error('No MC data in input argument')
end
C2.R=shiftdim(super.std(C2.RMC,0,1),1);
C2.RMC=[];
C2.MC=0;
end

function C2=var(C1)
%STD returns the standard deviation of the Monte Carlo data
C2=C1;
if isempty(C1.RMC)
   error('No MC data in input argument')
end
C2.R=shiftdim(super.std(C2.RMC,0,1),1).^2;
C2.RMC=[];
C2.MC=0;
end

function plot(varargin)
%PLOT plots the covariance function in covfun
%   plot(c1,c2,...,Property1,Value1,...)
%   Illustrates one or more covariance functions at the same time.
%
%   For covariance functions computed from estimated models, a confidence
%   bound can be added, based on Monte Carlo simulations of the covariance
%   function. The Monte Carlo data can also be shown in a scatter plot,
%   where all simulated random covariance functions are shown with half
%   the line width along with the nominal one.
%   The input objects m can be either LTI, SIG or COVFUN objects, but the first
%   one has to be a COVFUN object..
%   For SIG objects, sig2covfun is invoked, and similarly is lti2covfun
%   for LTI objects.
%
%   Property   Value        Description
%   ---------------------------------------------------------
%   MC        {20}          Monte Carlo simulations used in lti2covfun
%   conf      [{0},100]     Confidence level (0 mean no levels plotted)
%                           from MC data
%   scatter   'on'|{'off'}  Scatter plot of MC data
%   taumax    {30}          Maximum lag for which the covariance function
%                           is computed
%   interval  {'pos'}|'sym' Positive tau=0:taumax or
%                           symmetric tau=-taumax:taumax lag interval
%   axis      {gca}         Axis handle where plot is added
%   col       {'bgrmyk'}    Colors in order of appearance
%   fontsize  14            Font size
%   linewidth 2             Line width
%   Xlim      {}            Limits on x axis
%   Ylim      {}            Limits on y axis
%   legend    {}            Legend text
%
%   Examples:
%     m0=rand(arx(4));
%     y=simulate(m0,100);
%     mhat=estimate(arx(4),y);
%     csighat=covfun(y);
%     cmhat=covfun(mhat);
%     c0=covfun(m0);
%     plot(c0,csighat,cmhat)
%
%   See also:

%   Fredrik Gustafsson 08-May-2006
%   Copyright (c) 1994-2006 by COMSOL AB
%   $ Revision: 28-Oct-2019 $

N=0; K=0; optvar='';
k=0;
while k<length(varargin)
    k=k+1;
    if isstr(varargin{k})
        K=K+1;
        optvar{K}=varargin{k};
        K=K+1;
        k=k+1;
        optvar{K}=varargin{k};
    else
        N=N+1;
        M{N}=varargin{k};
    end
end

opt=struct('MC',20,'axis',0,'conf',90,'scatter','off','interval','pos','taumax',30,'col','bgrmyk','fontsize',14,'linewidth',2,'Xlim',[],'Ylim',[],'legend',[]);
opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold=ishold(opt.axis);

legon=0;
for k=1:N;
    kk=1+ceil(rem(k-1,length(opt.col)));
    m=M{k};
    if isa(m,'lss')
        c=ss2covfun(m,'taumax',opt.taumax,'MC',opt.MC);
        c=covfun(m);
    elseif isa(m,'covfun')
        c=m; % ok
    elseif isa(m,'sig')
        c=covfun(m); % ok
    else
        error('Objects must be either SS, SIG or COVFUN')
    end
    ny(k)=size(c);
    if k>1 & ny(k)~=ny(k-1),
        error('All COVFUN objects to plot must be of same size'),
    end
    R=c.R; tau=c.tau; RMC=c.RMC;
    if ~isempty(c.tlabel)
        ctlabel=c.tlabel;
    else
        ctlabel='tau';
    end
    if ~isempty(c.name)
        ctitle=c.name;
    else
        ctitle='';
    end
    if ~isempty(c.ylabel)
        cylabel=c.ylabel;
    else
        cylabel='R(tau)';
    end
    if strcmp(opt.interval,'sym');
        R=[R(end:-1:2,:,:);R]; tau=[-tau(end:-1:2);tau]; RMC=[RMC(end:-1:2,:,:,:);RMC];
    end
    for i=1:ny(k)
      for j=1:ny(k)
        if ny(k)>1,
            opt.axis=subplot(ny(k),ny(k),(i-1)*ny+j);
            oldhold = ishold(opt.axis);
        end
        if ~isempty(RMC)
            if strcmp(opt.scatter,'on')
                plot(tau,squeeze(RMC(:,:,i,j))',['-',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth/2);
                hold(opt.axis,'on')
            end
            if opt.conf>0 & opt.conf<=100
                p=confband(tau,squeeze(RMC(:,:,i,j))','plot',opt.col(kk),opt.axis,opt.conf);
                hold(opt.axis,'on')
                %if legon==0; legend([p1 p(1) p(3)],'Nominal','Confidence bound','Median','parent',opt.axis); legon=1; end
            end
        end
        p1=plot(tau,R(:,i,j),opt.col(kk),'parent',opt.axis,'linewidth',opt.linewidth);
        if k==N  % last plot
            if oldhold==1
                hold(opt.axis,'on')
            else
                hold(opt.axis,'off')
            end
            set(opt.axis,'xlabel',ctlabel,'ylabel',cylabel,'title',ctitle,'fontsize',opt.fontsize);
            if ny(k)>1
                set(opt.axis,'title',[ctitle,' E(y',num2str(i),'(t)y',num2str(j),'(t-tau)'])
            end
            if ~isempty(opt.Xlim)
                set(opt.axis,'Xlim',opt.Xlim);
            end
            if ~isempty(opt.Ylim)
                set(opt.axis,'Ylim',opt.Ylim);
            end
            if ~isempty(opt.legend)
                legend(opt.axis,opt.legend{:})
            end
        else
            hold(opt.axis,'on')
        end
      end
    end
end
end % plot
end % methods
end % class
