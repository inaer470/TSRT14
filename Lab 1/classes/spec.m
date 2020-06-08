classdef spec
%SPEC is the spectrum object
%   Fields:
%     Phi   is the spectrum. Dimension (nf,ny,ny)
%     f     is the frequency vector. Dimension nf
%     PhiMC contains MC realization. Dimension (MC,nf,ny,ny)
%
%   Examples:
%   t=(0:4095)';
%   f1=0.20; f2=0.24;
%   y=sin(2*pi*f1*t)+sin(2*pi*f2*t)+ndist(0,0.1);
%   M=30;
%   Phib=spec(y,'M',M,'method','blackman');
%   Phiw=spec(y,'M',M,'method','welch');
%   Phis=spec(y,'M',M,'method','smoothing');
%   semilogy(Phib,Phiw,Phis);
%   title('semilogy(Phib,Phiw,Phis)')

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $



properties (SetAccess = public)
    Phi; f; PhiMC;
    flabel; ylabel;
    fs; MC;
    name='empty'; desc; marker; markerlabel;
end

methods

function S=spec(varargin)
%SPEC constructor
%   S=spec               % empty object
%   S=spec(Phi,f)        % Phi(f,1:ny,1:ny)
%   S=spec(Phi,f,PhiMC)  % PhiMC(1:MC,f,1:ny,1:ny)
%   S=spec(s)           % Conversions from SIG and LTI (ARX, SS) objects
%   Optional public fields: fs MC name ylabel desc
%   Change these with R.fs=fs etc.
%
%   Example:
%     MC=30;
%     H=tf(1,[1 1.2 0.8],1);
%     Hf=freq(H);
%     f=Hf.f;
%     Phi=abs(Hf.H).^2;
%     PhiMC=4*(-0.5+rand(MC,size(Phi,1)))+repmat(Phi',MC,1);
%     Phi=spec(Phi,f,PhiMC);
%     plot(Phi);
%
%   See also: spec.plot, sig.estimate, lss.ss2spec, armax.armax2spec

%#
%#

Phi=[]; f=NaN; PhiMC=[];
fs=NaN; MC=0;
name='empty';  ylabel=[]; desc=[];
tlabel=''; marker=[]; markerlabel='';

if nargin==0,
    return % allow empty objects
elseif isa(varargin{1},'spec')  % Copy contents
    c=varargin{1};
elseif isa(varargin{1},'lss')
    c=ss2spec(varargin{1},varargin{2:end});
elseif isa(varargin{1},'arx')
    c=arx2spec(varargin{1},varargin{2:end});
elseif isa(varargin{1},'sig')
    c=sig2spec(varargin{1},varargin{2:end});
elseif isnumeric(varargin{1})
    Phi=varargin{1};
    if ~isnumeric(Phi), error('Phi must be numeric'), end
    if size(Phi,2)~=size(Phi,3)
        error('Phi must have dimension (nf,ny,ny)')
    end
    if nargin>1
        f=varargin{2};
        if ~isnumeric(f), error('Frequency vector f must be numeric'), end
        if length(f)~=size(Phi,1),
            error('Length of f must match the size of the first dimension of Phi')
        end
    end
    if nargin>2
        PhiMC=varargin{3};
        MC=size(PhiMC,1);
        if ~isempty(PhiMC) & (size(PhiMC,2)~=size(Phi,1) | size(PhiMC,3)~=size(Phi,2) | size(PhiMC,4)~=size(Phi,3))
            error('Last three dimensions of PhiMC must have size as Phi')
        end
    else
        PhiMC=[];
    end
else
   error('Incorrect syntax to SPEC')
end

if nargin>0 & (isa(varargin{1},'lti') | isa(varargin{1},'sig'))
%    Phi=c.Phi; f=c.f; PhiMC=c.PhiMC;
%    c=inherit(c,varargin{1},[' ',upper(class(varargin{1})),' -> SPEC']);
%    name=c.name; desc=c.desc; fs=c.fs; MC=c.MC;
%    ylabel=c.ylabel;
     S=c;
else
    S.Phi=Phi; S.f=f; S.PhiMC=PhiMC;
    S.fs=NaN; S.MC=0;
    S.name='Untitled';  S.ylabel=[]; S.desc=[];
end
end



function disp(S)
disp(['SPEC object: ',S.name])
end

function str=symbolic(S)
str=['SPEC: ',S.name];
end

function S=estimate(S,s,varargin)
%ESTIMATE estimates a spectrum from data
%   The sig method sig2spec is used
%   Example:
%   See also:
S=spec(s,varargin{:});
end

%==============================================
%-----OPERATORS FOR STOCHASTIC OBJECTS----------
%==============================================

function Sr=rand(S,n)
%RAND returns a cell array of SPEC realizations taken from MC data
%   Sr=rand(S,n);
%   n is the number of realizations (default n=1)
%   Sr is a cell array of length n

if nargin<2; n=1; end
if isempty(S.PhiMC)
   error('No MC data in input argument')
end
MC=size(S.PhiMC,1);
for k=1:n
   ind=ceil(MC*rand(1,1));
   Sr{k}=spec(shiftdim(S.PhiMC(ind,:,:),1),S.f);
   Sr{k}=inherit(Sr{k},S);
end

end

function S2=fix(S1)
%FIX removes the Monte Carlo data from the SPEC object
S2=S1;
S2.PhiMC=[];
S2.MC=0;
end

function S2=mean(S1)
%MEAN returns the mean signal of the Monte Carlo data
S2=S1;
if isempty(S1.PhiMC)
   error('No MC data in input argument')
end
S2.Phi=shiftdim(super.mean(S2.PhiMC,1),1);
S2.PhiMC=[];
S2.MC=0;
end

function S2=E(S1)
%E returns the mean signal of the Monte Carlo data
S2=mean(S1);
end

function S2=std(S1)
%STD returns the standard deviation of the Monte Carlo data
S2=S1;
if isempty(S1.PhiMC)
   error('No MC data in input argument')
end
S2.Phi=shiftdim(super.std(S2.PhiMC,0,1),1);
S2.PhiMC=[];
S2.MC=0;
end

function S2=var(S1)
%STD returns the standard deviation of the Monte Carlo data
S2=S1;
if isempty(S1.PhiMC)
   error('No MC data in input argument')
end
S2.Phi=shiftdim(super.var(S2.PhiMC,0,1),1);
S2.PhiMC=[];
S2.MC=0;
end


%==============================
%-----PLOTS--------------------
%==============================

function varargout=plot(varargin)
%PLOT calls specplot
%   plot(S1,S2,...,Property1,Value1,...)
%   See specplot for options
h=specplot(varargin{:});
if nargout > 0
  varargout{1} = h;
end
end

function varargout=semilogy(varargin)
%SEMILOGY calls specplot
%   semilogy(S1,S2,...,Property1,Value1,...)
%   See specplot for options
h=specplot(varargin{:},'plottype','semilogy');
if nargout > 0
  varargout{1} = h;
end
end

function varargout=semilogx(varargin)
%SEMILOGX calls specplot
%   semilogx(S1,S2,...,Property1,Value1,...)
%   See specplot for options
h=specplot(varargin{:},'plottype','semilogx');
if nargout > 0
  varargout{1} = h;
end
end

function varargout=loglog(varargin)
%LOGLOG calls specplot
%   loglog(S1,S2,...,Property1,Value1,...)
%   See specplot for options
h=specplot(varargin{:},'plottype','loglog');
if nargout > 0
  varargout{1} = h;
end
end

function varargout=specplot(varargin)
%SPECPLOT plots spectrum
%   specplot(spec1,spec2,...,Property1,Value1,...)
%   Illustrates one or more spectra at the same time.
%   spec can be the output from sig2spec, lti2spec
%   spec can also be a signal or LTI object, in which case
%   sig2spec and lti2spec are invoked with their default parameters.
%
%   Property  Value/{Default}      Description
%   ---------------------------------------------------------------------------
%   plottype  {'plot'} | 'semilogx'  Type of plot
%             'semilogy' | 'loglog'
%   MC        {30}                 Number of Monte Carlo simulations
%                                  used in lti2freq
%   conf      [{0},100]            Confidence level (0 means no levels plotted)
%                                  from MC data
%   scatter   'on' | {'off'}       Scatter plot of MC data
%   Xlim      [fmin fmax]          Focus on frequency axis
%   axis      {gca}                Axis handle where plot is added
%   linewidth {}                   Line width (default from global SIGNAL)
%   fontsize  {}                   Font size  (default from global SIGNAL)
%   conftype  {1} | 2              Confidence area (1) or lines (2)
%   legend    {'on'} | 'off'       Display spectrum data as legend
%   col       {'bgrmyck'}          Colors in order of appearance in spec1,spec2,...
%
%   Examples:
%     m0=randlti('arma',5);
%     y=lti2sig(m0,1000);
%     Phihat=sig2spec(y,'M',30,'method','welch');
%     mhat=sig2lti(y,'arma',5);
%     specplot(m0,Phihat,mhat);
%
%   See also: sig2spec, lti2spec

%   Fredrik Gustafsson 08-May-2006
%   Copyright (c) 1994-2006 by COMSOL AB
%   $ Revision: 28-Oct-2019 $


N=0; K=0; optvar='';
for k=1:length(varargin)
     if isstruct(varargin{k}) | isa(varargin{k},'spec')
        N=N+1;
        M{N}=varargin{k};
    else
        K=K+1;
        optvar{K}=varargin{k};
    end
end
opt=struct('bode','amplitude','plottype','plot','MC',30,'axis',0,'conf',90,'scatter','off','Xlim',[],'taumax',30,'linewidth',[],'fontsize',[],'conftype',1,'col','bgrmyk','legend','on');
opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

for k=1:N;
    kk=1+ceil(rem(k-1,length(opt.col)));
    m=M{k};
    if isa(m,'spec')
        % do nothing
    elseif isa(m,'lti') | isa(m,'sig')
        % convert
        m=spec(m);
    end
    tmp=symbolic(m);
    leg{k}=tmp(7:end);
    ny(k)=size(m.Phi,2);
    if k>1 & ny(k)~=ny(k-1)
        error('ny must be the same for all spectra')
    end
    for mm=1:ny(k)
        for nn=1:mm
            if ny(k)>1
                opt.axis=subplot(ny(k),ny(k),(mm-1)*ny(k)+nn); % magnitude
                oldhold = ishold(opt.axis);
                handle(k)=feval(opt.plottype,m.f,abs(m.Phi(:,mm,nn)),opt.col(kk),'parent',opt.axis,'linewidth',opt.linewidth);
                hold(opt.axis,'on')
                set(opt.axis,'xlabel','Frequency','ylabel','Power','fontsize',opt.fontsize);
                if nn==mm,
                    set(opt.axis,'title',['S',num2str(mm),num2str(nn)])
                else
                    set(opt.axis,'title',['|S',num2str(mm),num2str(nn),'|'])
                end
                if opt.conf>0 & opt.conf<=100 & isfield(m,'PhiMC') & ~isempty(m.PhiMC)
                    p=confband(m.f,abs(m.PhiMC(:,:,mm,nn)),opt.plottype,opt.col(kk),opt.axis,opt.conf,opt.conftype);
                end
                if strcmp(opt.scatter,'on') & isfield(m,'PhiMC')
                    feval(opt.plottype,m.f,abs(m.PhiMC(:,:,mm,nn)'),['-',opt.col(kk)],'parent',opt.axis,'linewidth',0.5);
                end
                if nn<mm  %phase
                    opt.axis=subplot(ny(k),ny(k),(nn-1)*ny(k)+mm); % magnitude
                    feval(opt.plottype,m.f,180/pi*angle(m.Phi(:,mm,nn)),opt.col(kk),'parent',opt.axis,'linewidth',opt.linewidth);
                    hold(opt.axis,'on')
                    set(opt.axis,'title',['arg(S',num2str(mm),num2str(nn),')'])
                    set(opt.axis,'xlabel','Frequency','ylabel','Power','fontsize',opt.fontsize);
                    if opt.conf>0 & opt.conf<=100 & isfield(m,'PhiMC') & ~isempty(m.PhiMC)
                        p=confband(m.f,180/pi*angle(m.PhiMC(:,:,mm,nn).'),opt.plottype,opt.col(kk),opt.axis,opt.conf,opt.conftype);
                    end
                    if strcmp(opt.scatter,'on') & isfield(m,'PhiMC')
                        feval(opt.plottype,m.f,180/pi*angle(m.PhiMC(:,:,mm,nn)'),['-',opt.col(kk)],'parent',opt.axis,'linewidth',0.5);
                    end
                    if k==N,
                       if oldhold==1
                           hold(opt.axis,'on')
                       else
                           hold(opt.axis,'off')
                       end
                    end
                end
            else  % Scalar signal, the standard case
                handle(k)=feval(opt.plottype,m.f,m.Phi,opt.col(kk),'parent',opt.axis,'linewidth',opt.linewidth);
                hold(opt.axis,'on')
                if opt.conf>0 & opt.conf<=100  & ~isempty(m.PhiMC)
                    p=confband(m.f,squeeze(m.PhiMC(:,:,1,1)'),opt.plottype,opt.col(kk),opt.axis,opt.conf,opt.conftype);
                end
                if  strcmp(opt.scatter,'on') & ~isempty(m.PhiMC)
                    feval(opt.plottype,m.f,squeeze(m.PhiMC(:,:,1,1)'),['-',opt.col(kk)],'parent',opt.axis,'linewidth',0.5);
                end
                xlabel('Frequency');
                ylabel('Power')
            end
            if k==N,
               hold(opt.axis,'off');
               if  strcmp(opt.legend,'on') & ny(k)==1
                  legend(handle,leg{:})
               end
            end
        end
    end
end
if ~isempty(opt.Xlim),
      set(opt.axis,'Xlim',opt.Xlim);
end

%legend('Nominal','Upper bound','Lower bound','Median')
if nargout>0
   varargout{1}=handle;
end
end

end % methods
end % class
