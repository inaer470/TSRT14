classdef ltimod
%LTIMOD is the parent of LSS, LTF, ARX, COV, SPEC, FREQ

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
% fs MC name xlabel ulabel ylabel desc
end

methods
function m=ltimod
end

function bode(varargin)
%BODE plots the Bode diagram of amplitude and phase
%   bode(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'view','bode','bodetype','both')
end

function bodeamp(varargin)
%BODEAMP plots the Bode diagram of amplitude only
%   bodeamp(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'view','bode','bodetype','amplitude')
end

function bodephase(varargin)
%BODEPHASE plots the Bode diagram of phase only
%   bodephase(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'view','bode','bodetype','phase')
end

function nyquist(varargin)
%NYQUIST plots the Nyquist curve
%   nyquist(s1,s2,...,Property1,Value1,...)
%   Note: confidence bounds not implemented for nyquist, only scatter plots
%   See help ltimod.plot for options
plot(varargin{:},'view','nyquist')
end

function zpplot(varargin)
%ZPPLOT plots the zeros and poles
%   zpplot(s1,s2,...,Property1,Value1,...)
%   Note: confidence bounds not implemented for zpplot, only scatter plots
%   See help ltimod.plot for options
plot(varargin{:},'view','zp')
end

function rlplot(varargin)
%RLPLOT plots the root locus
%   rlplot(s1,s2,...,Property1,Value1,...)
%   Note: Neither scatter plots nor confidence bounds implemented for rlplot
%   See help ltimod.plot for options
plot(varargin{:},'view','rl')
end

function loglog(varargin)
%LOGLOG plots the frequency response in logarithmic frequency and amplitude
%   loglog(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'plottype','loglog')
end

function semilogy(varargin)
%SEMILOGY plots the frequency response in logarithmic amplitude
%   semilogy(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'plottype','semilogy')
end

function semilogx(varargin)
%SEMILOGX plots the frequency response in logarithmic frequency
%   semilogx(s1,s2,...,Property1,Value1,...)
%   See help ltimod.plot for options
plot(varargin{:},'plottype','semilogx')
end

function plot(varargin)
%PLOT plots frequency response of LTIMOD objects
%   plot(obj1,obj2,...,Property1,Value1,...)
%   Illustrates various properties of LTIMOD (LSS/LTF) models and/or frequency models
%   This function over-loaded all classes belonging to the LTIMOD class
%
%   The following list of property-value pairs are available:
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   view      {'bode'}    Frequency response as Bode diagram
%             'nyquist'   Frequency response as Nyquist plot
%             'rootlocus' Roots of a+K*b with K as parameter
%             'zp'        Zero-pole plot
%             'all'       All of the above in four subplots
%   bodetype  {'amplitude'} Only amplitude response (default),
%             'both'      Both amplitude and phase response
%             'phase'     Only phase response
%   fmax      {'auto'}    Maximum frequency in bode and nyquist plots.
%   Kmax      {'auto'}    Maximum K in rlplot (K=[0:Kmax/1000:Kmax])
%   Kgrid     {'auto'}    Grid for K in rlplot as integer indeces
%   plottype  {'plot'} | 'semilogx' | 'semilogy' | 'loglog'
%   MC        {30}        Number of Monte Carlo simulations used in ltimod2freq
%   conf      [0,100] {0} Confidence level (default 0, no levels plotted) from MC data
%   conftype  {1} | 2     1=shaded confidence region, 2=dashed bounds and median
%   scatter   'on' | {'off'} Scatter plot of MC data
%   axis      {gca}       Axis handle where plot is added (does not apply
%                         for subplots in Bode diagrams of type 'both' or for MIMO)
%   col       {'bgrmyk'}  Colors in order of appearance in ltimod1,ltimod2,...
%   Xlim                  Limits on x axis
%   Ylim                  Limits on y axis
%   linewidth {}          Line width (default in global SIGNAL)
%   fontsize  {}          Font size  (default in global SIGNAL)
%   title     {'on'}|'off'Display a title of the plot
%   format    {'%11.2g'}  For text displays (K in rlplot, f in nyquist)

%
%   Example:
%      G1=randltimod('tf',3)
%      G2=exltimod('tf2d');
%      plot(G1,G2)
%   See also: ltimod, exltimod, randltimod


N=0; K=0; optvar='';
for k=1:length(varargin)
    if isa(varargin{k},'lss') |isa(varargin{k},'ltf') | isa(varargin{k},'ft') | isa(varargin{k},'arx')
        N=N+1;
        M{N}=varargin{k};
    else
        K=K+1;
        optvar{K}=varargin{k};
    end
end
opt=struct('view','bode','bodetype','amplitude','plottype','semilogy','MC',30,'axis',0,'conf',90,'Kmax','auto','fmax','auto','Kgrid','auto',...
    'scatter','off','Xlim',[],'Ylim',[],'taumax',30,'col','bgrmyk','Nellipse',20,'linewidth',[],'fontsize',[],'conftype',1,'title','on','format','%11.2f');

opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
% cla(opt.axis);
oldhold = ishold(opt.axis);


Ncol=length(opt.col);
if N>Ncol  % repeat the color pattern
   opt.col=char(kron(ones(1,ceil(N/Ncol)),opt.col));
end

abstype=lower(opt.plottype);
if strcmp(abstype,'plot')
    argtype='plot';
elseif strcmp(abstype,'semilogx')
    argtype='semilogx';
elseif strcmp(abstype,'semilogy')
    argtype='plot';
elseif strcmp(abstype,'loglog')
    argtype='semilogx';
else
    error('Incorrect plot type')
end

if strncmpi(opt.view,'bode',4)

for k=1:N;
  G=M{k};
  if isa(G,'ltf')
    G=freqfun(G,'MC',opt.MC,'fmax',opt.fmax);
  elseif isa(G,'lss')
    G=freqfun(G,'MC',opt.MC,'fmax',opt.fmax);
  elseif isa(G,'arx')
    G=freqfun(ltf(G),'MC',opt.MC,'fmax',opt.fmax);
  elseif isa(G,'ft')
    %ok
  else
    error('Input arguments must be either LSS, LTF, ARX or FT objects')
  end
  H=G.H; f=G.f; HMC=G.HMC;
  [nf,ny(k),nu(k)]=size(H);
  if k>1 % check dimensions
    if ny(k)~=ny(k-1)
      error('Objects have different output dimension')
    end
    if nu(k)~=nu(k-1)
      error('Objects have different input dimension')
    end
  end
  if strncmpi(opt.bodetype,'b',1)
    np=2;
  else
    np=1;
  end
  if ~isempty(G.ulabel)
      ulabel=G.ulabel;
  else  %Default signal names
      for i=1:nu(k)
          ulabel{i}=['u',num2str(i)];
      end
  end
  if ~isempty(G.ylabel)
      ylabel=G.ylabel;
  else  %Default signal names
      for i=1:ny(k)
          ylabel{i}=['y',num2str(i)];
      end
  end
  for j=1:ny(k)
    for i=1:nu(k)
      % Bode amplitude plot
      if strncmpi(opt.bodetype,'a',1) | strncmpi(opt.bodetype,'b',1)
        if nu(k)>1 | ny(k)>1 | np>1
            opt.axis=subplot(np*ny(k),nu(k),(np*(j-1))*nu(k)+i);
            oldhold = ishold(opt.axis);
        end
        feval(abstype,f,abs(H(:,j,i)),opt.col(k),'parent',opt.axis,'linewidth',opt.linewidth)
        hold(opt.axis,'on')
        set(opt.axis,'xlabel','Frequency [Hz]','ylabel','Magnitude','fontsize',opt.fontsize)
        if ~isempty(HMC)
          if strcmp(opt.scatter,'on')
            feval(abstype,f,abs(squeeze(HMC(:,:,j,i)))',['-',opt.col(k)],'parent',opt.axis,'linewidth',0.5);
          end
          if opt.conf>0 & opt.conf<=100
            p=confband(f,abs(squeeze(HMC(:,:,j,i)))',abstype,opt.col(k),opt.axis,opt.conf,opt.conftype);
          end
        end
        if (nu(k)>1 | ny(k)>1) | (~isempty(G.ulabel) | ~isempty(G.ylabel))
            title([ulabel{i},' to ',ylabel{j}])
        end
        if k==N
           if oldhold==1
              hold(opt.axis,'on')
           else
              hold(opt.axis,'off')
           end
        end
        if ~isempty(opt.Xlim),
           set(opt.axis,'Xlim',opt.Xlim)
        end
        if ~isempty(opt.Ylim),
           set(opt.axis,'Ylim',opt.Ylim)
        end
    end
    % Bode phase
    if strncmpi(opt.bodetype,'p',1) | strncmpi(opt.bodetype,'b',1)
        if nu(k)>1 | ny(k)>1 | np>1
            opt.axis=subplot(np*ny(k),nu(k),(np*j-1)*nu(k)+i);
            oldhold = ishold(opt.axis);
%        else
%            opt.axis=subplot(np*ny(k),nu(k),(np*j-1)*nu(k)+i);
        end
        feval(argtype,f,180/pi*angle(H(:,j,i)),opt.col(k),'parent',opt.axis,'linewidth',opt.linewidth)
        set(opt.axis,'xlabel','Frequency [Hz]','ylabel','Phase [deg]','fontsize',opt.fontsize)
        hold(opt.axis,'on')
        if ~isempty(HMC)
            if strcmp(opt.scatter,'on')
                feval(argtype,f,180/pi*angle(squeeze(HMC(:,:,j,i)))',['-',opt.col(k)],'parent',opt.axis,'linewidth',0.5);
            end
            if opt.conf>0 & opt.conf<=100
                p=confband(f,180/pi*angle(squeeze(HMC(:,:,j,i)))',argtype,opt.col(k),opt.axis,opt.conf,opt.conftype);
            end
        end
        if (nu(k)>1 | ny(k)>1) | (~isempty(G.ulabel) | ~isempty(G.ylabel))
            title([ulabel{i},' to ',ylabel{j}])
        end
      end
      if k==N
         if oldhold==1
            hold(opt.axis,'on')
         else
            hold(opt.axis,'off')
         end
      end
      if ~isempty(opt.Xlim),
        set(opt.axis,'Xlim',opt.Xlim)
      end
      if ~isempty(opt.Ylim),
         set(opt.axis,'Ylim',opt.Ylim)
      end
  end
end
end
if strcmp('opt.title','on')
    set(opt.axis,'title','Bode')
end


elseif strncmpi(opt.view,'ny',2)  % Nyquist curve
    for k=1:N;
        G=M{k};
        if isa(G,'ltf')
            G=freqfun(G,'MC',opt.MC,'fmax',opt.fmax);
        elseif isa(G,'lss')
            G=freqfun(G,'MC',opt.MC,'fmax',opt.fmax);
        elseif isa(G,'arx')
            G=freqfun(ltf(G),'MC',opt.MC,'fmax',opt.fmax);
        elseif isa(G,'ft')
            %ok
        else
            error('Input arguments must be either LSS, LTF, ARX or FT objects')
        end
        H=G.H; f=G.f; HMC=G.HMC;
        [nf,ny(k),nu(k)]=size(H);
        if k>1 % check dimensions
          if ny(k)~=ny(k-1)
            error('Objects have different output dimension')
          end
          if nu(k)~=nu(k-1)
            error('Objects have different input dimension')
          end
        end
        if ~isempty(G.ulabel)
            ulabel=G.ulabel;
        else  %Default signal names
            for i=1:nu(k)
                ulabel{i}=['u',num2str(i)];
            end
        end
        if ~isempty(G.ylabel)
            ylabel=G.ylabel;
        else  %Default signal names
            for i=1:ny(k)
                ylabel{i}=['y',num2str(i)];
            end
        end

        for j=1:ny(k)
          for i=1:nu(k)
              if nu(k)>1 | ny(k)>1
                  opt.axis=subplot(ny(k),nu(k),(j-1)*nu(k)+i);
                  oldhold = ishold(opt.axis);
              end
              plot(real(H(:,j,i)),imag(H(:,j,i)),opt.col(k),'parent',opt.axis,'linewidth',opt.linewidth)
              text(real(H(1,j,i)),imag(H(1,j,i)),['f=',num2str(f(1),opt.format)],'parent',opt.axis,'fontsize',opt.fontsize)
              text(real(H(end,j,i)),imag(H(end,j,i)),['f=',num2str(f(end),opt.format)],'parent',opt.axis,'fontsize',opt.fontsize)
              hold(opt.axis,'on')
              if ~isempty(HMC)
                  if strcmp(opt.scatter,'on')
                      plot(real(squeeze(HMC(:,:,j,i)))',imag(squeeze(HMC(:,:,j,i)))',['-',opt.col(k)],'parent',opt.axis,'linewidth',0.5)
                  end
                  if opt.conf>0 & opt.conf<=100
                      ind=round(linspace(1,length(f),opt.Nellipse));
      %                confellipse(real([H(ind,j,i) squeeze(HMC(:,ind,j,i))]),imag([H(ind,j,i) squeeze(HMC(:,ind,j,i))]),opt.col(k),opt.axis,opt.conf);
                  end
              end
              set(opt.axis,'xlabel','Real','ylabel','Imag','fontsize',opt.fontsize)
              if (nu(k)>1 | ny(k)>1) | (~isempty(G.ulabel) | ~isempty(G.ylabel))
                  title([ulabel{i},' to ',ylabel{j}])
              end
              if strcmp('opt.title','on')
                  set(opt.axis,'title','Nyquist')
              end
              if ~isempty(opt.Xlim),
                  set(opt.axis,'Xlim',opt.Xlim)
              end
              if ~isempty(opt.Ylim),
                  set(opt.axis,'Ylim',opt.Ylim)
              end
              phi=linspace(0,2*pi,20);
              drawnow
              ax=axis;
              plot(-1+(ax(2)-ax(1))/50*sin(phi),(ax(4)-ax(3))/50*cos(phi),'r','parent',opt.axis,'linewidth',opt.linewidth);
              drawnow
              ax=axis;
              plot(ax(1:2),[0 0],'k','parent',opt.axis,'linewidth',opt.linewidth)
              plot([0 0],ax(3:4),'k','parent',opt.axis,'linewidth',opt.linewidth)
              if k==N
                 if oldhold==1
                     hold(opt.axis,'on')
                 else
                     hold(opt.axis,'off')
                 end
              end
          end
      end
    end
elseif strncmpi(opt.view,'zp',2)    % Zero-pole
    for k=1:N;
        G=M{k};
        if isa(G,'ltf')
            [z,p,K,zmc,pmc,Kmc]=zpk(G);
        elseif isa(G,'arx')
            [z,p,K,zmc,pmc,Kmc]=zpk(ltf(G));
        elseif isa(G,'lss')
            [z,p,K,zmc,pmc,Kmc]=zpk(G);
        else
            error('Input arguments must be either LSS, LTF or ARX objects')
        end
        ny(k)=size(z,1);
        nu(k)=size(z,3);
        if k>1 % check dimensions
          if ny(k)~=ny(k-1)
            error('Objects have different output dimension')
          end
          if nu(k)~=nu(k-1)
            error('Objects have different input dimension')
          end
        end
        if ~isempty(G.ulabel)
            ulabel=G.ulabel;
        else  %Default signal names
            for i=1:nu(k)
                ulabel{i}=['u',num2str(i)];
            end
        end
        if ~isempty(G.ylabel)
            ylabel=G.ylabel;
        else  %Default signal names
            for i=1:ny(k)
                ylabel{i}=['y',num2str(i)];
            end
        end

        for j=1:ny(k)
          for i=1:nu(k)
              if nu(k)>1 | ny(k)>1
                  opt.axis=subplot(ny(k),nu(k),(j-1)*nu(k)+i);
                  oldhold = ishold(opt.axis);
              end
              plot(real(p),imag(p),['*',opt.col(k)],'parent',opt.axis,'linewidth',opt.linewidth);
              hold(opt.axis,'on')
              plot(real(z(j,:,i)),imag(z(j,:,i)),['o',opt.col(k)],'parent',opt.axis,'linewidth',opt.linewidth);
              if strcmp(opt.scatter,'on') & ~isempty(pmc)
                  if ~isempty(pmc)
                    plot(real(pmc(:,:)'),imag(pmc(:,:)'),['.',opt.col(k)],'parent',opt.axis,'linewidth',0.5);
                  end
                  hold(opt.axis,'on')
                  if ~isempty(zmc)
                     plot(real(squeeze(zmc(:,j,:,i))),imag(squeeze(zmc(:,j,:,i))),['.',opt.col(k)],'parent',opt.axis,'linewidth',0.5);
                  end
              end
              axis('equal');
              phi=linspace(0,2*pi,200);
              if  isempty(G.fs) | G.fs>0 % fix
                  plot(cos(phi),sin(phi),'k','parent',opt.axis,'linewidth',opt.linewidth) % Unit circle
              end
              drawnow
              a=axis;
              plot(a(1:2),[0 0],'k','parent',opt.axis,'linewidth',opt.linewidth)
              plot([0 0],a(3:4),'k','parent',opt.axis,'linewidth',opt.linewidth)
              drawnow
              axis(a);
              if k==N
                 if oldhold==1
                     hold(opt.axis,'on')
                 else
                     hold(opt.axis,'off')
                 end
              end
              set(opt.axis,'xlabel','Real','ylabel','Imag','fontsize',opt.fontsize)
              if ~isempty(opt.Xlim),
                  set(opt.axis,'Xlim',opt.Xlim)
              end
              if ~isempty(opt.Ylim),
                  set(opt.axis,'Ylim',opt.Ylim)
              end
              if (nu(k)>1 | ny(k)>1) | (~isempty(G.ulabel) | ~isempty(G.ylabel))
                  title([ulabel{i},' to ',ylabel{j}])
              end
          end
        end
    end

elseif strncmpi(opt.view,'rl',2)  % root locus
    for k=1:N;
        G=M{k};
        if ~( isa(G,'ltf') | isa(G,'lss')| isa(G,'arx'))
            error('RLPLOT: Input arguments must be either LSS, LTF  or ARX objects')
        end
        if strcmp(opt.Kmax,'auto')
            Gtmp=ltf(G);
            if isnan(Gtmp.fs)
              indb=find(Gtmp.b(1,:,1)~=0);
              Kmax=abs(10/Gtmp.b(1,indb(1),1));
            else
              Kmax=10^(ceil(log10(abs(10*sum(Gtmp.a)/sum(Gtmp.b(1,:,1))))));
            end
        else
            Kmax=opt.Kmax;
            if ~isnumericscalar(Kmax) | Kmax<=0
                error('RLPLOT: Kmax must be a positive scalar')
            end
        end
        K=0:Kmax/1000:Kmax;
        G=lss(G);
        z=rl(G,K);  % Use only SS form which works for MIMO
        NK=length(K);
        if strcmp(opt.Kgrid,'auto')
            Kgrid=1:ceil((NK-1)/1):NK;
        else
            Kgrid=opt.Kgrid;
        end
    %axis([-2 2 -2 2])
        for j=1:size(z,1)
            plot(real(z(j,:)),imag(z(j,:)),['-',opt.col(k)],'parent',opt.axis,'linewidth',opt.linewidth);
            hold(opt.axis,'on')
            for l=1:length(Kgrid)
                plot(real(z(:,Kgrid(l))),imag(z(:,Kgrid(l))),['o',opt.col(k)],'parent',opt.axis,'linewidth',opt.linewidth);
                text(real(z(:,Kgrid(l))),imag(z(:,Kgrid(l))),['K=',num2str(K(Kgrid(l)),opt.format)],'color',opt.col(k),'parent',opt.axis);
            end
        end
        axis('equal');
        drawnow
        ax=axis;
        plot(ax(1:2),[0 0],'k','parent',opt.axis,'linewidth',opt.linewidth/2)
        plot([0 0],ax(3:4),'k','parent',opt.axis,'linewidth',opt.linewidth/2)
        if G.fs>0
            phi=linspace(0,2*pi,200);
            plot(cos(phi),sin(phi),'k','parent',opt.axis,'linewidth',opt.linewidth/2) % Unit circle
        end
    end
    set(opt.axis,'xlabel','Real','ylabel','Imag','fontsize',opt.fontsize)
    if strcmp('opt.title','on')
        set(opt.axis,'title','Root locus')
    end
    if oldhold==1
        hold(opt.axis,'on')
    else
        hold(opt.axis,'off')
    end
    drawnow
    if ~isempty(opt.Xlim),
        set(opt.axis,'equal','off','Xlim',opt.Xlim)
    end
    if ~isempty(opt.Ylim),
        set(opt.axis,'equal','off','Ylim',opt.Ylim)
    end
elseif strncmpi(opt.view,'all',3)
     subplot(2,2,1), plot(M{:},optvar{:},'view','bode')
     subplot(2,2,2), plot(M{:},optvar{:},'view','nyquist')
     subplot(2,2,3), plot(M{:},optvar{:},'view','zp')
     subplot(2,2,4), plot(M{:},optvar{:},'view','root')
end
end %plot
end %methods
end %ltimod
