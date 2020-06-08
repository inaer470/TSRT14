classdef ft

%FT Fourier Transform class for signals
%   An object in the FT class represents a signal by its Fourier Transform
%   (Y,f). It also represents uncertainties in the transform with MC samples.
%
%   Fields    Description
%   ft.Y      The mandatory primary data field
%   ft.f      Mandatory frequency
%   ft.YMC    Monte Carlo realizations
%   ft.MC     Number of MC data
%   ft.ylabel Labels for the output
%   ft.flabel Labels for the frequency
%   ft.marker Discrete events that are marked in the signal
%   ft.markerlabel Labels for the events
%
%   help ft.ft    gives help for the constructor
%   methods ft    lists all methods for this class
%
%   Example:
%   t=(0:31)';                   % Time interval
%   f1=0.22;                     % Frequency
%   y=sin(2*pi*f1*t);            % Sinusoid signal
%   Y=fft(y);                    % Y=DFT(y)
%   Y1=ft(Y(1:16),(0:15)/32);    % Direct definition of FT object
%   Y2=ft(sig(y,t));             % Conversion from SIG to FT object
%   plot(Y1,Y2)

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
    Y; f; YMC;
    flabel; ylabel;
    MC;
    name; desc; marker; markerlabel;
end

methods
function Yf=ft(varargin)
%FT constructor
%   Y=ft(y)        conversion from a SIG, LTF or LSS object
%   Y=ft(Y,f,YMC)
%   Y   is the Fourier transform (approximation) of a signal
%       For multi-variate signals, Y is (nf,ny) dimensional
%   f   is the frequency vector in Y(f)
%   YMC contains Monte Carlo realizations of Y stored in a matrix
%       which is of size (MC,nf,ny)

%#
%#


ylabel=[]; flabel=[];
ulabel=[]; xlabel=[]; marker=[]; markerlabel=[];
name='';  desc=[];
MC=0; YMC=[]; U=[];

if nargin<1
    return % allow empty objects
end
if isa(varargin{1},'ft')  % Do nothing
    Yf=varargin{1};
    Y=Yf.Y; f=Yf.f; YMC=Yf.YMC; MC=Yf.MC;
    name=Yf.name; desc=Yf.desc;
    ylabel=Yf.ylabel; flabel=Yf.flabel;

elseif isa(varargin{1},'sig')  % Convert from SIG object
    z=varargin{1};
    Yf=sig2ft(z,varargin{2:end});
    Y=Yf.Y; f=Yf.f;  YMC=Yf.YMC; MC=Yf.MC;
    name=Yf.name; desc=Yf.desc;
    ylabel=Yf.ylabel; flabel=Yf.flabel;

elseif isa(varargin{1},'ltf')  % Convert from SIG object
    z=varargin{1};
    Yf=ltf2ft(z,varargin{2:end});
    Y=Yf.Y; f=Yf.f;  YMC=Yf.YMC; MC=Yf.MC;
    name=Yf.name; desc=Yf.desc;
    ylabel=Yf.ylabel; flabel=Yf.flabel;

elseif isa(varargin{1},'lss')  % Convert from SIG object
    z=varargin{1};
    Yf=ltf2ft(ltf(z),varargin{2:end});
    Y=Yf.Y; f=Yf.f;  YMC=Yf.YMC; MC=Yf.MC;
    name=Yf.name; desc=Yf.desc;
    ylabel=Yf.ylabel; flabel=Yf.flabel;

else
    if nargin<2
        error('Input to FT must either be a SIG object, or Y,f (and YMC)')
    end
    Y=varargin{1};
    f=varargin{2};
    if length(f)~=size(Y,1)
        error('Length of f and number of rows in Y must be the same')
    end
    if nargin>2;
        YMC=varargin{3};
        MC=size(YMC,1);
        if MC>0 & (size(Y,1)~=size(YMC,2) | size(Y,2)~=size(YMC,3))
            error('Last two dimensions of YMC must be the same as the dimensions of Y')
        end
    end
end
    Yf.Y=Y; Yf.f=f; Yf.YMC=YMC; Yf.MC=MC;
    Yf.name=name; Yf.desc=desc;
    Yf.ylabel=ylabel; Yf.flabel=flabel;
end



function out=arrayread(varargin)
%ARRAYREAD used to select frequencies and sub-signals
%   Y2=Y(f,yind);
%   f is a vector of frequency indices
%   yind is a vector of signal indices

nf=size(Y,1);
ny=size(Y,2);
freqind=varargin{1};
if isstr(freqind)
  if strcmp(freqind,':')
    freqind=1:nf;
  else
    error('Inputs to arrayread must be a vector or :')
  end
end
if any(freqind>length(f))
    error('Requested index larger than length of frequency vector')
end
if any(freqind<1)
    error('Frequency index must be positive')
end
if any(freqind~=round(freqind))
    error('Frequency index must be integers')
end
if nargin>1
  yind=varargin{2};
  if isstr(yind)
    if strcmp(yind,':')
        yind=1:ny;
    else
       error('Signal index must be a vector or :')
    end
  end
else
  yind=1:ny;
end
if any(yind>ny)
     error('Signal index larger than the number of signals')
end
if any(yind<1)
    error('Signal index must be positive integers')
end
if any(yind~=round(yind))
    error('Signal index must be integers')
end

YMC2=[];
if ~isempty(YMC)
    YMC2=YMC(:,freqind,yind);
end
out=ft(Y(freqind,yind),f(freqind),YMC2);
end

function [nf,ny,MC]=size(Yf);
%SIZE returns the sizes
%   [nf,ny,MC]=size(Yf);
nf=size(Yf.Y,1);
ny=size(Yf.Y,2);
MC=size(Yf.YMC,1);
end

function disp(Y)
%DISP
if isempty(Y.Y)
   disp('Empty FT object')
   return
end
disp(['FT object: ',Y.name])
disp(['  Number of frequency points: ',num2str(length(Y.f))])
disp(['  Number of outputs: ',num2str(size(Y.Y,2))])
disp(['  Number of Monte Carlo samples: ',num2str(Y.MC)])
end

%=============================
%-----Inverse-----------------
%=============================

function y=ift(Y)
%IFT is the inverse Fourier Transform
%   y=ift(Y)
%   The function basically computes the inverse FFT and then
%   truncates the result to compensate for zero padding.
%
%   Limitation: The FT must have been computed on a uniform frequency grid
%
%   Example:
%     yv=sin((1:100)/3);
%     y=sig(yv,3);
%     Y=ft(y);
%     yi=ift(Y)

if any(diff(diff(Y.f)~=0))
   error('FT.IFT: frequency grid must be uniform')
end
fmax=max(Y.f);
fs=2*fmax;

Yv=Y.Y;
[Nf,ny]=size(Yv);
Yv=[Yv;conj(Yv(end-1:-1:2,:))];
for k=1:ny
   y(:,k)=ifft(Yv(:,k));
   ind1=find(abs(y)<1e-13*max(abs(y)));
   if ~isempty(ind1)
      ind2=find(diff([0;ind1(:)])~=1);
      ind(k)=ind1(ind2(end))-1;
   else
      ind(k)=0;
   end
end
tmax=max(ind);
if tmax>0
    y=y(1:tmax,:);  % Truncate last part where all signals are zero
end
yMC=[];
if ~isempty(Y.YMC)
    for k=1:size(Y.MC,1)
        ytmp=ift(ft(shiftdim(Y.YMC(k,:,:),1),fs));
        yMC(k,:,:)=ytmp.y;
     end
end
y=sig(y,fs,yMC);
y=inherit(y,Y,'FT -> SIG');
end


%=============================
%-----OPERATIONS--------------
%=============================

function Y2=abs(Y1);
%ABS computes the absolute value of a FT object
%   Y2=abs(Y1);
Y2=ft(abs(Y1.Y),Y1.f,abs(Y1.YMC));
Y2.ylabel=Y1.ylabel;
end

function Y2=angle(Y1);
%ANGLE computes the angle (phase) of a FT object
%   Y2=angle(Y1);
Y2=ft(angle(Y1.Y),Y1.f,angle(Y1.YMC));
Y2.ylabel=Y1.ylabel;
end
function Y2=real(Y1);
%REAL computes the real value of a FT object
%   Y2=real(Y1);
Y2=ft(real(Y1.Y),Y1.f,real(Y1.YMC));
Y2.ylabel=Y1.ylabel;
end
function Y2=imag(Y1);
%IMAG computes the imaginary value of a FT object
%   Y2=imag(Y1);
Y2=ft(imag(Y1.Y),Y1.f,imag(Y1.YMC));
Y2.ylabel=Y1.ylabel;
end
function Y=plus(Y1,Y2);
%PLUS adds two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('plus',Y1,Y2);
end
function Y=minus(Y1,Y2);
%MINUS subtracts two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('minus',Y1,Y2);
end
function Y=times(Y1,Y2);
%TIMES multiplies two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('times',Y1,Y2);
end
function Y=mtimes(Y1,Y2);
%MTIMES multiplies two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('times',Y1,Y2);
end
function Y=rdivide(Y1,Y2);
%RDIVIDE divides two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('rdivide',Y1,Y2);
end
function Y=divide(Y1,Y2);
%DIVIDE divides two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('rdivide',Y1,Y2);
end
function Y=mrdivide(Y1,Y2);
%MRDIVIDE divides two Fourier transforms, or one FT and a scalar/vector
Y=evalfun2('rdivide',Y1,Y2);
end

function Y=evalfun2(op,arg1,arg2)
%EVALFUN2 contains common code for plus, minus, times and divide

if isa(arg1,'ft')
    Y1=arg1;
    Y2=arg2;
else
    Y1=arg2;
    Y2=arg1;
end
if isnumericscalar(Y2)  % op on a constant
    Y=Y1;
    Y.Y=feval(op,Y.Y,Y2);
    Y.YMC=feval(op,Y.YMC,Y2);
elseif isnumeric(Y2)    % op on a vector
    if size(Y1.Y)==size(Y2)
        Y=Y1;
        Y.Y=feval(op,Y.Y,Y2);
        if Y.MC>0 & ~isempty(Y.YMC);
            tmp(1,:,:)=Y2;
            Y2MC=repmat(tmp,[Y.MC 1 1]);
            Y.YMC=feval(op,Y.YMC,Y2MC);
        end
    else
        error(['Using ',upper(op),' on a vector/matrix and a FT object, size(Y.Y) must coincide with the size of the vector/matrix'])
   end
elseif isa(Y2,'ft') % op on two signals
    if size(Y1.Y)==size(Y2.Y) & size(Y1.f)==size(Y2.f)
        Y=Y1;
        Y.Y=feval(op,Y.Y,Y2.Y);
        if Y1.MC==Y2.MC
           % Nothing to do
        elseif Y1.MC>0 & Y2.MC>0 & ~isempty(Y1.YMC) & ~isempty(Y2.YMC)
           if Y1.MC~=Y2.MC
               MC=min([Y1.MC,Y2.MC]);
               Y.YMC=feval(op,Y1.YMC(1:MC,:,:),Y2.YMC(1:MC,:,:));
           else
               Y.YMC=feval(op,Y1.YMC,Y2.YMC);
           end
        elseif Y1.MC==0 & ~isempty(Y2.YMC);
            MC=Y2.MC;
            tmp(1,:,:)=Y1.Y;
            Y1MC=repmat(tmp,[MC 1 1]);
            Y.YMC=feval(op,Y2.YMC,Y1MC);
        elseif Y2.MC==0 & ~isempty(Y1.YMC);
            MC=Y1.MC;
            tmp(1,:,:)=Y2.Y;
            Y2MC=repmat(tmp,[MC 1 1]);
            Y.YMC=feval(op,Y1.YMC,Y2MC);
        else
            error('FT: Check that the MC number coincides with the actual number of MC simulations in the YMC field')
        end
    else
        error(['Using ',upper(op),'on two FT objects, N1==N2, fs1==fs2, f1==f2'])
   end
else
   error(['Inappropriate argument for ',upper(op),' in FT objects'])
end
end




%==============================================
%-----OPERATORS FOR STOCHASTIC OBJECTS----------
%==============================================
function Yfr=rand(Yf,n)
%RAND returns a cell array of FT realizations taken from MC data
%   Yfr=rand(Yf,n);
%   n is the number of realizations (default n=1)
%   Yfr is a cell array of length n

if nargin<2; n=1; end
if isempty(Yf.YMC)
   error('No MC data in input argument')
end
MC=size(Yf.YMC,1);
for k=1:n
   ind=ceil(MC*rand(1,1));
   Yfr{k}=ft(shiftdim(Yf.YMC(ind,:,:),1),Yf.f);
   Yfr{k}=inherit(Yfr{k},Yf);
end

end
function Yf2=fix(Yf1)
%FIX removes the Monte Carlo data from the FT object
Yf2=Yf1;
Yf2.YMC=[];
Yf2.MC=0;
end
function Yf2=mean(Yf1)
%MEAN returns the mean signal of the Monte Carlo data
Yf2=Yf1;
if isempty(Yf1.YMC)
   error('No MC data in input argument')
end
Yf2.Y=shiftdim(super.mean(Yf2.YMC,1),1);
Yf2.YMC=[];
Yf2.MC=0;
end
function Yf2=E(Yf1)
%E returns the mean signal of the Monte Carlo data
Yf2=mean(Yf1);
end
function Yf2=std(Yf1)
%STD returns the standard deviation of the Monte Carlo data
Yf2=Yf1;
if isempty(Yf1.YMC)
   error('No MC data in input argument')
end
Yf2.Y=shiftdim(super.std(Yf2.YMC,0,1),1);
Yf2.YMC=[];
Yf2.MC=0;
end

function Yf2=var(Yf1)
%STD returns the standard deviation of the Monte Carlo data
Yf2=Yf1;
if isempty(Yf1.YMC)
   error('No MC data in input argument')
end
Yf2.Y=shiftdim(super.std(Yf2.YMC,0,1),1).^2;
Yf2.YMC=[];
Yf2.MC=0;
end




%=============================
%-----PLOT FUNCTIONS----------
%=============================
function plot(varargin)
%PLOT illustrates Fourier Transforms in a plot
%   plot(Y1,Y2,...,Property1,Value1,...);
%   This function is the plot method for FT objects
%   Multiple plots of several FT's are possible.
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   type       1,{2},3     Interval for f,
%                          1:[0,fs], 2:[0,fs/2], 3:[-fs/2,fs/2]
%   plottype  'plot'|{'semilogy'}|'semilogx'|'loglog'
%   conf      {0}         Confidence level in per cent
%   conftype  1|2          Band or limits
%   scatter   'on'|{'off'} Scatter plot of Monte Carlo data
%   Xlim      {}           Limits on x axis
%   Ylim      {}           Limits on y axis
%   axis      {gca}        Axis handle where plot is added
%   col       {'bgrmyk'}   Colors in order of appearance
%   fontsize  {14}         Font size
%   linewidth {2}          Line width
%
%   Example:
%     load genera
%     yd=detrend(y1,3);
%     Y=ft(yd);
%     plot(Y)   % Default view
%     plot(Y,'Xlim',[-0.02 0.02],'plottype','plot')
%
%   See also: freq.plot

%   Fredrik Gustafsson 08-May-2006
%   Copyright (c) 1994-2006 by COMSOL AB
%   $ Revision: 28-Oct-2019 $

N=0; K=0; optvar='';
k=0;
while k<length(varargin)
    k=k+1;
    if isa(varargin{k},'ft')
        N=N+1;
        M{N}=varargin{k};
    elseif isstr(varargin{k})
            %Property value pair
            K=K+1;
            optvar{K}=varargin{k};
            K=K+1;
            k=k+1;
            optvar{K}=varargin{k};
    else   % Data vector
            error('Inputs must be a FT object, or property value pairs')
    end
end

opt=struct('type',3,'axis',0,'col','bgrmyk','fontsize',14,'linewidth',2,'Xlim',[],'Ylim',[],'plottype','plot','conftype',1,'conf',0,'scatter','off');
opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

for k=1:N;
    Yf=M{k};
    if ~isreal(Yf.Y), Yf=abs(Yf); end  % absolute value is default
    [Ny,ny(k)]=size(Yf);
    Y=Yf.Y;
    f=Yf.f(:);
    YMC=Yf.YMC;
    if k>1 % check dimensions
      if ny(k)~=ny(k-1)
        error('Signals have different output dimension')
      end
    end
    if ~isempty(Yf.ylabel)
       if ny(k)==1 & isstr(Yf.ylabel),
          ylabel{1}=Yf.ylabel;
       else
          ylabel=Yf.ylabel;
       end
    else  %Default signal names
       for i=1:ny(k)
          ylabel{i}=['y',num2str(i)];
       end
    end
    if ~isempty(Yf.flabel)
       flabel=Yf.flabel;
    else
       flabel='Frequency';
    end
    switch opt.type
            case 1
                f=[f(1:Ny); f(Ny)+f(2:Ny-1)];
                Y=[Y(1:Ny,:); Y(Ny-1:-1:2,:)];
                if ~isempty(YMC)
                    YMC=cat(2,YMC(:,1:Ny,:), YMC(:,Ny-1:-1:2,:));
                end
            case 2
                f=f(1:Ny);
                Y=Y(1:Ny,:);
                if ~isempty(YMC)
                    YMC=YMC(:,1:Ny,:);
                end
            case 3
                f=[f(2:Ny-1)-f(Ny); f(1:Ny)];
                Y=[Y(Ny-1:-1:2,:); Y(1:Ny,:)];
                if ~isempty(YMC)
                    YMC=cat(2,YMC(:,Ny-1:-1:2,:), YMC(:,1:Ny,:));
                end
    end
    for j=1:ny(k)
        if ny(k)>1
          opt.axis=subplot(ny(k),1,j);
          oldhold = ishold(opt.axis);
        end
        eval([opt.plottype,'(f,abs(Y(:,j)),opt.col(k),''parent'',opt.axis,''linewidth'',opt.linewidth);'])
        hold(opt.axis,'on')
        if ~isempty(YMC)
           if strcmp(opt.scatter,'on')  % scatter always interpolated
              plot(f,YMC(:,:,j)',['--',opt.col(k)],'parent',opt.axis,'linewidth',0.5*opt.linewidth);
           end
           if opt.conf>0 & opt.conf<=100
              p=confband(f,YMC(:,:,j)','plot',opt.col(k),opt.axis,opt.conf);
           end
        end
        xlabel(opt.axis, flabel);
        if ny(k)==1 & ~isempty(Yf.name)
           set(opt.axis,'title',Yf.name);
        end

        if k==N  % last signal
         drawnow
         ylim=get(opt.axis,'Ylim');
         set(opt.axis,'Ylim',[ylim(1)-0.1*diff(ylim) ylim(2)+0.1*diff(ylim)])
         if ~isempty(opt.Xlim),
            set(opt.axis,'Xlim',opt.Xlim)
         end
         if ~isempty(opt.Ylim),
            set(opt.axis,'Ylim',opt.Ylim)
         end
         title(ylabel{j})
         if oldhold==1
           hold(opt.axis,'on')
         else
           hold(opt.axis,'off')
         end
        end
    end
end
set(opt.axis,'fontsize',opt.fontsize);
end

function semilogy(varargin)
%SEMILOGY calls ft.plot
%   semilogy(Y1,Y2,...,Property1,Value1,...)
%   See ft.plot for options
plot(varargin{:},'plottype','semilogy');
end
function semilogx(varargin)
%SEMILOGX calls ft.plot
%   semilogx(Y1,Y2,...,Property1,Value1,...)
%   See ft.plot for options
plot(varargin{:},'plottype','semilogx');
end

function loglog(varargin)
%LOGLOG calls ft.plot
%   loglog(Y1,Y2,...,Property1,Value1,...)
%   See ft.plot for options
plot(varargin{:},'plottype','loglog');
end %loglog

end % methods
end % class
