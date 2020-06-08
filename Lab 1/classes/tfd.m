classdef tfd
%TFD is the class for time-frequency descriptions of signals
%
%   help tfd.tfd    gives help for the constructor
%   methods tfd     lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
Y, t, f
fs, MC, name, xlabel, ulabel, ylabel, tlabel, desc, marker, markerlabel
end

methods

function Yt=tfd(varargin)
%TFD is the Time-Frequency Description object with constructions:
%   Yt=tfd          Empty object
%   Yt=tfd(Y,t,f)   Direct construction
%   Yt=tfd(mt)      Conversion from RARX model
%   Yt=tfd(z)       Estimation from signal object, equivalent to
%                   Yt=estimate(tfd,z) and Yt=sig2tfd(z)
%
%   The matrix in the field Y must have dimension (length(f),length(t))

Y=[]; t=[]; f=[]; fs=[]; MC=[]; name=[]; fs=1;
xlabel=[]; ulabel=[]; ylabel=[]; tlabel=[]; desc=[]; marker=[]; markerlabel=[];
if nargin==0  % empty object
elseif isa(varargin{1},'sig')
    Yt=sig2tfd(varargin{1},varargin{2:end});
    Y=Yt.Y; t=Yt.t; f=Yt.f;
    name=Yt.name; xlabel=Yt.xlabel; ulabel=Yt.ulabel; ylabel=Yt.ylabel;
    desc=Yt.desc; marker=Yt.marker; markerlabel=Yt.markerlabel;
elseif isa(varargin{1},'rarx')
    Yt=rarx2tfd(varargin{1},varargin{2:end});
    Y=Yt.Y; t=Yt.t; f=Yt.f;
    name=Yt.name; xlabel=Yt.xlabel; ulabel=Yt.ulabel; ylabel=Yt.ylabel;
    desc=Yt.desc; marker=Yt.marker; markerlabel=Yt.markerlabel;
elseif isa(varargin{1},'tfd')
    Yt=varargin{1};
    Y=Yt.Y; t=Yt.t; f=Yt.f;
    name=Yt.name; xlabel=Yt.xlabel; ulabel=Yt.ulabel; ylabel=Yt.ylabel;
    desc=Yt.desc; marker=Yt.marker; markerlabel=Yt.markerlabel;
elseif nargin==3
    Y=varargin{1};
    t=varargin{2};
    f=varargin{3};
    if size(Y,1)~=length(f)
        error('TF constructor: number of rows in E must equal the length of f')
    end
    if size(Y,2)~=length(t)
        error('TF constructor: number of columns in E must equal the length of t')
    end
    fs=2; MC=0; name='';
    xlabel=[]; ulabel=[]; ylabel=[]; tlabel=[]; desc=[]; marker=[]; markerlabel=[];
else
    error('Incorrect syntax to TFD constructor')
end
Yt.Y=Y; Yt.t=t; Yt.f=f; Yt.fs=fs; Yt.MC=MC; Yt.name=name;
Yt.xlabel=xlabel; Yt.ulabel=ulabel; Yt.ylabel=ylabel; Yt.tlabel=tlabel;
Yt.desc=desc; Yt.marker=marker; Yt.markerlabel=markerlabel;

end % constructor


function out=fieldread(arg)
out=eval(arg);
end

function fieldwrite(arg1,arg2)
error('Cannot change protected fields in TFD object.')
end

function disp(Yt)
disp(['TFD object: ',Yt.name])
end

function Yt=estimate(Yt,s,varargin)
%ESTIMATE estimates a TFD from data
%   The sig method sig2tfd is used
%   Example:
%   See also:
Yt=tfd(s,varargin{:});
end

function plot(Yt,varargin)
%PLOT illustrates a time frequency descriptions (TFD) as an image
%   plot(Yt,Property1,Value1,...)
%   See tfdplot for options
tfdplot(Yt,varargin{:},'view','image')
end

function image(Yt,varargin)
%IMAGE illustrates a time frequency descriptions (TFD) as an image
%   image(Yt,Property1,Value1,...)
%   See tfdplot for options
tfdplot(Yt,varargin{:},'view','image')
end

function surf(Yt,varargin)
%SURF illustrates a time frequency descriptions (TFD) as a surf plot
%   surf(Yt,Property1,Value1,...)
%   See tfdplot for options
tfdplot(Yt,varargin{:},'view','surf')
end

function contour(Yt,varargin)
%CONTOUR illustrates a time frequency descriptions (TFD) as a contour plot
%   contour(Yt,Property1,Value1,...)
%   See tfdplot for options
tfdplot(Yt,varargin{:},'view','contour')
end

function mesh(Yt,varargin)
%MESH illustrates a time frequency descriptions (TFD) as a mesh plot
%   mesh(Yt,Property1,Value1,...)
%   See tfdplot for options
tfdplot(Yt,varargin{:},'view','mesh')
end

function tfdplot(Yt,varargin)
%TFDPLOT illustrates a time frequency descriptions (TFD)
%   tfdplot(Yt,Property1,Value1,...)
%   Different 3D views of the TFD object are posssible.
%   Note that only one TFD object at the time can be shown.
%   The function is called from ltv2tfd, sig2tfd and ltvplot.
%
%   Yt       The TFD object
%
%   Property  Value/{default}  Description
%   axis      {gca}            Axis handle where plot is added
%   view      {'contour'}      TFD as contour plot
%              'surf'          TFD as surf plot
%              'mesh'          TFD as mesh plot
%              'image'         TFD as an image
%   histeq    {'on'} | 'off'   Histogram equalization for energy (z) values
%   linewidth {}               Line width (default from global SIGNAL)
%   fontsize  {}               Font size  (default from global SIGNAL)
%
%   Examples:
%     m=exrarx('ar2');
%     Yt=tfd(m);
%     surf(Yt)
%
%   See also: tfd, sig.sig2tfd,

%   Fredrik Gustafsson 08-May-2006
%   Copyright (c) 1994-2006 by COMSOL AB
%   $ Revision: 28-Oct-2019 $

opt=struct('axis',0,'linewidth',[],'fontsize',[],'view','image','histeq','on');
opt=optset(opt,varargin);
if opt.axis==0; opt.axis=gca; end

Y=Yt.Y; t=Yt.t; f=Yt.f;
if strcmpi(opt.histeq,'on')
    Y=histeq(Y);  % Histogram equalization
end
[tt,ff]=meshgrid(t,f);
cla(opt.axis)
if strcmpi(opt.view,'surf')
    surf(t,f,Y,'parent',opt.axis)
elseif strcmpi(opt.view,'mesh')
    mesh(t,f,Y,'parent',opt.axis)
elseif strcmpi(opt.view,'image')
    imagesc([t(1) t(end)],[f(1) f(end)],Y(end:-1:1,:),'parent',opt.axis)
    set(opt.axis,'Ydir','normal')
else
  contour(tt,ff,log(Y),'parent',opt.axis,'linewidth',opt.linewidth)
end
if ~isempty(Yt.name)
   title(Yt.name)
end
set(opt.axis,'fontsize',opt.fontsize);
xlabel('Time')
ylabel('Frequency')
zlabel('abs(H(t,f))')
set(opt.axis,'box','on')
end

end %methods
end %tfd
