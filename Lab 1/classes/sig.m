classdef sig
%SIG is the data class for signals
%   An object in the SIG class stores numerical data (y,x,u) from a system,
%   where y is the mandatory measurement, u is an optional input and
%   x is an optional state. It also stores uncertainties represented with
%   covariances or samples. Other fields are sampling frequency and/or time
%   vector, and names for all variables.
%
%   Fields     Description
%   sig.y      The mandatory primary data field
%   sig.fs     Mandatory sampling frequency, fs=NaN for non-uniform data
%   sig.u      Optional input that generated the output
%   sig.x      Optional state in the system that generated the output
%   sig.Py     Covariance for output
%   sig.Px     Covariance for state
%   sig.yMC    Monte Carlo realizations of the output
%   sig.xMC    Monte Carlo realizations of the state
%   sig.MC     Number of MC data
%   sig.ylabel Labels for the output
%   sig.ulabel Labels for the input
%   sig.xlabel Labels for the state
%   sig.tlabel Labels for the time
%   sig.marker Discrete events that are marked in the signal
%   sig.markerlabel Labels for the events
%
%   help sig.sig    gives help for the constructor
%   methods sig     lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
    y; u; x; t; yMC; xMC; fs;
    tlabel; xlabel; ulabel; ylabel; Px; Py;
    nn; MC;
    name; desc; marker; markerlabel; moviefile; userdata;
end

methods
function s=sig(varargin)
%SIG constructor
%   sig(y,fs)    % Uniformly sampled time series y[k]=y(k/fs)
%   sig(y,t)     % Non-uniformly sampled time series y(t) representing cont time
%   sig(y,t,u)   % Input to io system
%   sig(y,t,u,x) % State vector in state space system
%   Covariances of state and output can be provided as
%   sig(y,fs,u,x,Py,Px)
%   MC data can be provided as a matrix if MC > ny and MC > nx:
%   sig(y,fs,u,x,yMC,xMC)
%
%   Sizes of the arguments
%   signal  dim1 dim2 dim3
%   t         N
%   y         N   ny
%   u         N   nu
%   x         N   nx
%   yMC       MC   N   ny
%   xMC       MC   N   nx
%   Py        N   ny   ny
%   Px        N   nx   nx
%
%   Generally, the second input argument denotes sampling frequency is scalar
%   and time if a vector. When defining one sample, the second argument is time.
%
%   MC controls the number of Monte Carlo simulations, while the other ones
%   are used in displays and plots.
%
%   protected fields: y u x t yMC xMC fs  xlabel ulabel ylabel Px Py
%   public fields: MC name tlabel desc marker markerlabel moviefile userdata
%
%   Examples:
%     fs=1; N=32; MC=30;
%     t1=(0:1:N-1)'*fs;
%     y1=sin(0.3*t1);
%     z1=sig(y1,fs)     % Constructor for uniformly sampled signals
%     t2=sort(N*rand(31,1));
%     y2=sin(0.3*t2);
%     z2=sig(y2,t2)     % Constructor for non-uniformly sampled signals
%     plot(z1), hold on
%     stem(z2), hold off
%     z1n=z1+0.1*ndist(0,1);                  % Implicit definition of MC data
%     yMC=repmat(y1',MC,1)+0.1*randn(MC,N);   % Explicit definition
%     z1n=sig(y1,fs,[],[],yMC);               % of Monte Carlo data
%     plot(z1n,'conf',90)

%#
%#


MC=30;
if nargin<1
    return % allow empty objects
    %error('SIG constructor requires at least one input argument.')
end
y=varargin{1};
if isa(y,'struct')
   try
     s.y=y.y;
   catch
     error('SIG constructor: the input is a struct that does not have the required fields y')
   end
   try s.u=y.u;  catch u=[]; end
   try s.x=y.x;  catch x=[]; end
   s.fs=1;
   try s.fs=y.fs;  end
   try s.t=y.t;  end
   if isempty(s.t)
       N=size(y.y,1);
       s.t=(0:N-1)'/s.fs;
   end
   try s.yMC=y.yMC;  end
   try s.xMC=y.xMC;  end
   try s.xlabel=y.xlabel;  end
   try s.ulabel=y.ulabel;  end
   try s.ylabel=y.ylabel;  end
   try s.Px=y.Px;  end
   try s.Py=y.Px; end
   try s.MC=y.MC; end
   try s.name=y.name;  end
   try s.tlabel=y.tlabel;  end
   try s.desc=y.desc; end
   try s.marker=y.marker;  end
   try s.markerlabel=y.markerlabel;  end
   try  s.moviefile=y.moviefile;  end
   try s.userdata=y.userdata; end
   s.nn=[size(s.x,2) size(s.u,2), size(s.y,2)];
   return
end

if isa(y,'sig')
    %error('SIG constructor input is already a SIG object')
   s.y=y.y; s.u=y.u; s.x=y.x; s.t=y.t; s.yMC=y.yMC; s.xMC=y.xMC; s.fs=y.fs;
   s.xlabel=y.xlabel; s.ulabel=y.ulabel; s.ylabel=y.ylabel;
   s.Px=y.Px; s.Py=y.Px;
   s.MC=y.MC;
   s.name=y.name; s.tlabel=y.tlabel; s.desc=y.desc;
   s.marker=y.marker; s.markerlabel=y.markerlabel;
   s.moviefile=y.moviefile; s.userdata=y.userdata;
   s.nn=y.nn;
   return
end

if isa(y,'nl') % NL model
   x=y.x0';
   P=cov(y.px0);
   s.y=x; s.x=x;
   s.Px=P; s.Py=P;
   s.u=[]; s.t=1; s.yMC=[]; s.xMC=[]; s.fs=y.fs;
   s.nn=[size(s.x,2) 0, size(s.y,2)];
   s.xlabel=y.xlabel; s.ulabel=[]; s.ylabel=y.xlabel;
   s.MC=0;
   s.name=y.name; s.tlabel=y.tlabel; s.desc=y.desc;
   %   s.nn=y.nn;
   return
end

if any(isnan(y)) | any(isinf(y))
%    error('SIG: Data y contains NaN or Inf')
end
if  ~isnumeric(y) % | isempty(y) allow empty y from 061214
    error('SIG: y must be a numeric vector or matrix.')
end
if size(y,1)==1 & nargin<3;
%   y=y.';  % User friendly fix too dangerous
end
[N,ny]=size(y);
t=(0:N-1)';
fs=1;
if nargin>1
    if isnumericscalar(varargin{2})
        if N>1
	   fs=varargin{2};
           if fs<=0
               error('SIG constructor: fs must be positive')
           end
           t=(0:N-1)'/fs;
        else
           t=varargin{2};
        end
    elseif isvector(varargin{2});
        t=varargin{2};
        if length(t)~=size(y,1)
            error('SIG constructor: Length of t must equal number of rows in y')
        end
        fs=NaN;
    else
        error('SIG constructor: Incorrect second input to SIG (t or fs expected)')
    end
else
    fs=1;
    t=(0:N-1)';
end
if nargin>2
    u=varargin{3};
    if isempty(u)
       u=zeros(N,0);
    end
    if size(u,1)~=N
       error('SIG constructor: size(u,1) must equal size(y,1)=N')
    end
else
       u=zeros(N,0);
end
if nargin>3
    x=varargin{4};
    if isempty(x)
       x=zeros(N,0);
    end
    if size(x,1)~=N
       error('SIG constructor: size(x,1) must equal size(y,1)=N')
    end
else
    x=zeros(N,0);
end
if nargin>4 & (size(varargin{5},2)~=ny  | size(varargin{5},1)~=N )
    yMC=varargin{5};
    MC=size(yMC,1);
    if ~isempty(yMC) & size(yMC,2)~=N
       error('SIG constructor: size(yMC,2) must equal size(y,1)=N')
    end
    if ~isempty(yMC) & size(yMC,3)~=size(y,2)
       error('SIG constructor: size(yMC,3) must equal size(y,2)=ny')
    end
    Py=[];
elseif nargin>4 & size(varargin{5},2)==ny
    Py=varargin{5};
    if ~isempty(Py)
       if size(Py,1)~=N
          error('SIG constructor: size(Py,1) must equal size(y,1)=N')
       end
       if size(Py,3)~=size(Py,2)
          error('SIG constructor: Py(i,:,:) must be square')
       end
       if ~iscov(squeeze(Py(1,:,:)))
          error('SIG constructor: Py(i,:,:) must be covariance matrices')
       end
    end
    yMC=[];
else
    yMC=[];
    Py=[];
end
if nargin>5 & size(varargin{6},2)~=size(varargin{6},3)
    xMC=varargin{6};
    if ~isempty(xMC) & size(xMC,1)~=MC
        error('SIG constructor: The number of Monte Carlo samples in xMC and yMC should be the same')
    end
    if ~isempty(xMC) & size(xMC,2)~=N
       error('SIG constructor: size(xMC,2) must equal size(y,1)=N')
    end
    if ~isempty(xMC) & size(xMC,3)~=size(x,2)
       error('SIG constructor: size(xMC,3) must equal size(x,2)=nx')
    end
else
    xMC=[];
end
if nargin>5 & size(varargin{6},2)==size(varargin{6},3)
    Px=varargin{6};
    if ~isempty(Px)
       if size(Px,1)~=N
          error('SIG constructor: size(Px,1) must equal size(x,1)=N')
       end
       if size(Px,3)~=size(Px,2)
          error('SIG constructor: Px(i,:,:) must be square')
       end
       if ~iscov(squeeze(Px(1,:,:)))
          error('SIG constructor: Px(i,:,:) must be covariance matrices')
       end
    end
else
    Px=[];
end
ylabel=[]; ulabel=[]; xlabel=[];
for i=1:size(y,2);
      ylabel{i}=['y',num2str(i)];
end
for i=1:size(u,2);
      ulabel{i}=['u',num2str(i)];
end
for i=1:size(x,2);
      xlabel{i}=['x',num2str(i)];
end
   s.y=y; s.u=u; s.x=x; s.yMC=yMC; s.xMC=xMC; s.t=t; s.fs=fs;
   s.nn=[size(x,2) size(u,2), size(y,2)];
   s.xlabel=xlabel; s.ulabel=ulabel; s.ylabel=ylabel;
   s.Px=Px; s.Py=Px;
   s.MC=MC;
end



%function y=get.y(s)
%  y=s.y;
%end


function s=set.y(s,y)
   if ~isnumeric(y)
      error('SIG: y must be numeric.')
   end
   if isempty(s.y) | isequal(size(s.y),size(y)) | isnan(s.fs)
      s.y=y;
   else
     error('SIG: Cannot change size of y by assignment')
   end
end
function s=set.u(s,u)
   if ~isnumeric(u)
      error('SIG: u must be numeric.')
   end
   if isempty(s.u) | isequal(size(s.u),size(u))
      s.u=u;
   else
     error('SIG: Cannot change size of u by assignment')
   end
end

function s=set.x(s,x)
   if ~isnumeric(x)
      error('SIG: x must be numeric.')
   end
   if isempty(s.x) | isequal(size(s.x),size(x))
      s.x=x;
   else
     error('SIG: Cannot change size of x by assignment')
   end
end

function s=set.yMC(s,yMC)
   if ~isnumeric(yMC)
      error('SIG: yMC must be numeric.')
   end
   if isempty(s.yMC) | isequal(size(s.yMC),size(yMC))
      s.yMC=yMC;
   else
     error('SIG: Cannot change size of yMC by assignment')
   end
end

function s=set.xMC(s,xMC)
   if ~isnumeric(xMC)
      error('SIG: xMC must be numeric.')
   end
   if isempty(s.xMC) | isequal(size(s.xMC),size(xMC))
      s.xMC=xMC;
   else
     error('SIG: Cannot change size of xMC by assignment')
   end
end

function s=set.t(s,t)
   if ~isnumeric(t)
      error('SIG: t must be numeric.')
   end
   if ~isreal(t)
      error('SIG: t must be real')
   end
   [n,m]=size(t);
   if n>1 & m>1
      error('SIG: t must be a vector')
   end
   t=t(:);
   if length(t)~=size(s.y,1)
      error('SIG: t must be of the same length as y')
   end
   if any(diff(t)<0)
%      error('SIG: t must be non-decreasing')
   end
   s.t=t;
   Ts=diff(t);
   if isempty(Ts) | any(abs(diff(Ts))>10*eps)
      s.fs=NaN;
   else
      s.fs=Ts(1);
   end
end

function s=set.fs(s,fs)
   if ~isnumeric(fs)
      error('SIG: fs must be numeric.')
   end
   if isnumericscalar(fs) & fs>0
      s.fs=fs;
%      disp('SIG warning: cannot change time indeces to match this fs, change t recommended')
   elseif isnumericscalar(fs) & fs==0
      s.fs=NaN;
   elseif isnan(fs)
      s.fs=NaN;
   else
      error('SIG: Sampling frequency must be positive')
   end
end

function s=set.xlabel(s,xlabel)
   if ~iscell(xlabel) & isstr(xlabel) & s.nn(1)==1
      s.xlabel={xlabel};   % Cellify
   elseif iscell(xlabel) & length(xlabel)==s.nn(1)
      s.xlabel=xlabel;
   elseif isempty(xlabel)
      for i=1:s.nn(1)
         s.xlabel{i}=['x',num2str(i)];
      end
   else
      error(['SIG fieldwrite: Field xlabel must be a cell of length nx, [], or string if nx=1'])
   end
end
function s=set.ylabel(s,ylabel)
   if ~iscell(ylabel) & isstr(ylabel) & s.nn(3)==1
      s.ylabel={ylabel};   % Cellify
   elseif iscell(ylabel) & length(ylabel)==s.nn(3)
      s.ylabel=ylabel;
   elseif isempty(ylabel)
      for i=1:s.nn(3)
         s.ylabel{i}=['y',num2str(i)];
      end
   else
      error(['SIG fieldwrite: Field ylabel must be a cell of length ny, [], or string if ny=1'])
   end
end

function s=set.ulabel(s,ulabel)
   if ~iscell(ulabel) & isstr(ulabel) & s.nn(2)==1
      s.ulabel={ulabel};   % Cellify
   elseif iscell(ulabel) & length(ulabel)==s.nn(2)
      s.ulabel=ulabel;
   elseif isempty(ulabel)
      for i=1:s.nn(2)
         s.ulabel{i}=['u',num2str(i)];
      end
   else
      error(['SIG fieldwrite: Field ulabel must be a cell of length nu, [], or string if nu=1'])
   end
end

function s=set.tlabel(s,tlabel)
   if isempty(tlabel) | isstr(tlabel)
      s.tlabel=tlabel;
   else
      error(['SIG fieldwrite: Field tlabel must be a string'])
   end
end

function s=set.name(s,name)
   if isempty(name) | isstr(name)
      s.name=name;
   else
      error(['SIG fieldwrite: Field name must be a string'])
   end
end

function s=set.desc(s,desc)
   if isempty(desc) |isstr(desc)
      s.desc=desc;
   else
      error(['SIG fieldwrite: Field desc must be a string'])
   end
end

function s=set.moviefile(s,moviefile)
   if isempty(moviefile) | isstr(moviefile)
      s.moviefile=moviefile;
   else
      error(['SIG fieldwrite: Field moviefile must be a string'])
   end
end

function s=set.markerlabel(s,markerlabel)
   if isempty(markerlabel) | isstr(markerlabel)
      s.markerlabel=markerlabel;
   else
      error(['SIG fieldwrite: Field markerlabel must be a string'])
   end
end

function s=set.marker(s,marker)
   if isempty(marker) | isvector(marker)
      s.marker=marker;
   else
      error(['SIG fieldwrite: Field marker must be a vector'])
   end
end


function ind = end(obj, k, n)
% END returns last index in the array
%  ind = end(obj, k, n)
%  obj the object
%  k   the index in the expression using the end syntax
%  n   the total number of indices in the expression
% ind  the index value to use in the expression

  ind = size(obj, k);
end


%function zout=index(s,S)
function zout=subsref(s,S)
%SUBSREF is used to pick out subsignals from SIG objects
%   stij=s(t,i,j)
%   t is a vector of time indices
%   i is a vector of index/indices, corresponding to the outputs
%   j is a vector of index/indices, corresponding to the inputs
%   Example: z(:,2,3) gives the vector signal from input 2 to output 3

%S.type,S.subs
if length(S)>1
        zout = s.(S(1).subs)(S(2:end).subs{:});
else
switch S.type
case '.'
     zout =s.(S.subs);
case '{}'
     error('SIG: sig arrays not defined')
case '()'
     tind=S.subs{1};
     if length(S.subs)<3
         j=':';
     else
       j=S.subs{3};
     end
     if length(S.subs)<2
         i=':';
     else
       i=S.subs{2};
     end
     [N,ny,nu,nx]=size(s);
     if isstr(tind)
         if strcmp(tind,':')
             tind=(1:N)';
         else
            error('Inputs to arrayread must be a vector or :')
         end
     end
     if isstr(i)
         if strcmp(i,':')
             i=1:ny;
         else
            error('Inputs to arrayread must be a vector or :')
         end
     end
     if isstr(j)
         if strcmp(j,':')
             j=1:nu;
         else
            error('Inputs to arrayread must be a vector or :')
         end
     end
     if any(i>ny)
          error('Output (row) index larger than the number of outputs')
     end
     if any(j>nu)
          error('Input (column) index larger than the number of inputs')
     end
     if any(tind<1) | any(tind>N)
         error('Specified time interval out of range')
     end
     yy=s.y(tind,i);
     tt=s.t(tind);
     if ~isempty(s.u)
         uu=s.u(tind,j);
     else
         uu=[];
     end
     if ~isempty(s.x)
         xx=s.x(tind,:);
     else
         xx=[];
     end
     if ~isempty(s.yMC)
         yyMC=s.yMC(:,tind,i);
     else
         yyMC=[];
     end
     if ~isempty(s.xMC)
         xxMC=s.xMC(:,tind,:);
     else
         xxMC=[];
     end
     if ~isempty(s.Py)
         Pyy=s.Py(tind,i,i);
     else
         Pyy=[];
     end
     if ~isempty(s.Px)
         Pxx=s.Px(tind,:,:);
     else
         Pxx=[];
     end
     zout=sig(yy,tt,uu,xx,yyMC,xxMC);
     zout.Py=Pyy;
     zout.Px=Pxx;
     zout.fs=s.fs;
     zout.name=s.name;
     zout.xlabel=s.xlabel;
     if ~isempty(s.ulabel)
        zout.ulabel={s.ulabel{j}};
     end
     if ~isempty(s.ylabel)
        zout.ylabel = {s.ylabel{i}};
     end
     zout.tlabel=s.tlabel;
     zout.desc=s.desc;
     zout.marker=s.marker;
     zout.markerlabel=s.markerlabel;
end  % switch
end  % if
end  % subsref

function list(s)
%LIST lists available benchmark signals with a brief description
%   list(sig)
pp=which('genera.mat');
help(pp(1:end-11))
end

function z=horzcat(z1,z2);
%HORZCAT concatenates two SIG objects to larger output dimension
%   z=horzcat(z1,z2) or z=[z1 z2]
%   The outputs are stacked together as [z1.y z2.y]
%   Requirements:
%      The outputs must have the same length
%      For input-output relations, the inputs must be the same
%        In the case of different inputs, use append
%      The sampling rates must be equal
%   The user properties are inherited from the first argument

if size(z1,1)~=size(z2,1)
   error('The outputs must be of the same length')
elseif ~isequal(z1.u,z2.u)
   error('The inputs must be the same')
elseif ~isequalwithequalnans(z1.fs,z2.fs)
   error('The sampling rates must be the same')
else
  z = struct('u', z1.u,...
             'x', [z1.x, z2.x],...
             'y', [z1.y, z2.y],...
             'fs', z1.fs,...
             't', z1.t,...
             'yMC', cat(3, z1.yMC, z2.yMC),...
             'xMC', cat(3, z1.xMC, z2.xMC),...
             'xlabel', z1.xlabel,...
             'ulabel', z1.ulabel,...
             'ylabel', z1.ylabel,...
             'Px', z1.Px,...
             'Py', z1.Py,...
             'MC', z1.MC,...
             'name', z1.name,...
             'tlabel', z1.tlabel,...
             'desc', z1.desc,...
             'marker', z1.marker,...
             'markerlabel', z1.markerlabel,...
             'moviefile', z1.moviefile,...
             'userdata', z1.userdata);
  z = sig(z);
end
end

function z=vertcat(z1,z2);
%VERTCAT concatenates two SIG objects in time
%   z=vertcat(z1,z2)  or  z=[z1;z2]
%   The outputs are stacked together as [z1.y z2.y]
%   Requirements:
%      The outputs and inputs must have the same dimensions
%      The sampling rates must be equal
%   The user properties are inherited from the first argument

if size(z1,3)~=size(z2,3)
   error('The inputs must be of the same dimension')
elseif size(z1,2)~=size(z2,2)
   error('The outputs must be of the same dimension')
elseif size(z1,4)~=size(z2,4)
   error('The states must be of the same dimension')
elseif ~isequalwithequalnans(z1.fs,z2.fs)
   error('The sampling rates must be the same')
else
  z = struct('u', [z1.u; z2.u],...
             'x', [z1.x; z2.x],...
             'y', [z1.y; z2.y],...
             'fs', z1.fs,...
             't', z1.t,...
             'yMC', cat(2, z1.yMC, z2.yMC),...
             'xMC', cat(2, z1.xMC, z2.xMC),...
             'xlabel', z1.xlabel,...
             'ulabel', z1.ulabel,...
             'ylabel', z1.ylabel,...
             'Px', z1.Px,...
             'Py', z1.Py,...
             'MC', z1.MC,...
             'name', z1.name,...
             'tlabel', z1.tlabel,...
             'desc', z1.desc,...
             'marker', z1.marker,...
             'markerlabel', z1.markerlabel,...
             'moviefile', z1.moviefile,...
             'userdata', z1.userdata);
   if isnan(z1.fs) % time vector
      z.t=[z1.t(:); z1.t(end)+z2.t(:)];
   end
   z = sig(z);
end
end

function z=append(z1,z2);
%APPEND concatenates two SIG objects to MIMO signals
%   z=append(z1,z2)
%   The outputs and inputs are stacked side by side
%   as [z1.y z2.y] and [z1.u z2.u]
%   Requirements:
%      The outputs must have the same length
%      For input-output relations, the inputs should be different
%        In the case of identical inputs, use vertcat
%      For time series, append is the same as vertcat
%      The sampling rates must be equal
%   The user properties are inherited from the first argument

if size(z1,1)~=size(z2,1)
   error('The outputs must be of the same length')
elseif isequal(z1.u,z2.u)
   error('The inputs must be different, otherwise use vertcat')
elseif ~isequalwithequalnans(z1.fs,z2.fs)
   error('The sampling rates must be the same')
else
  z = struct('u', [z1.u, z2.u],...
             'x', [z1.x, z2.x],...
             'y', [z1.y, z2.y],...
             'fs', z1.fs,...
             't', z1.t,...
             'yMC', z1.yMC,...
             'xMC', z1.xMC,...
             'xlabel', z1.xlabel,...
             'ulabel', z1.ulabel,...
             'ylabel', z1.ylabel,...
             'Px', z1.Px,...
             'Py', z1.Py,...
             'MC', z1.MC,...
             'name', z1.name,...
             'tlabel', z1.tlabel,...
             'desc', z1.desc,...
             'marker', z1.marker,...
             'markerlabel', z1.markerlabel,...
             'moviefile', z1.moviefile,...
             'userdata', z1.userdata);
   z = sig(z);
end
end

function disp(s)
%DISP function for SIG objects
if isempty(s.y)
   disp('Empty SIG object')
   return
end
y=s.y; u=s.u; x=s.x; nn=s.nn; fs=s.fs;
name=s.name; desc=s.desc;
MC=s.MC; yMC=s.yMC;
moviefile=s.moviefile;

N=size(y,1);
str2=['N = ',num2str(N)];
ny=nn(3);
str2=[str2,',  ny = ',num2str(ny)];
if isnan(fs)
   str3=['continuous time '];
else
   fsstr=['(fs = ',num2str(fs),') '];
   str3=['discrete time ',fsstr];
end
if isempty(u) & isempty(x)
   str='time series';
elseif ~isempty(u) & isempty(x)
   str='input-output data';
   nu=size(u,2);
   str2=[str2,', nu = ',num2str(nu)];
elseif ~isempty(u) & ~isempty(x)
   str='input-output state space data';
   nx=size(x,2);
   nu=size(u,2);
   str2=[str2,', nu = ',num2str(nu)];
   str2=[str2,', nx = ',num2str(nx)];
elseif isempty(u) & ~isempty(x)
   str='stochastic state space data (no input)';
   nx=size(x,2);
   str2=[str2,', nx = ',num2str(nx)];
end
disp(['SIG object with ',str3,str])
if ~isempty(name)
   disp(['  Name:        ',name])
end
if ~isempty(desc)
    [nd,md]=size(desc);
    bl=blanks(15*(nd-1));
    disp([ ['  Description: ';reshape(bl,nd-1,15)] desc])
end
disp(['  Sizes:       ',str2])
if MC>0
    disp(['  MC is set to: ',num2str(MC)])
    disp(['  #MC samples:  ',num2str(size(yMC,1))])
end
if ~isempty(moviefile)
    disp(['Movie file:  ',moviefile])
end
end

function info(z)
%INFO displays user specified information of a SIG object
disp(['Name:    ',z.name])
if ~isempty(z.desc)
    disp(['Description: '])
    disp(z.desc)
end
disp('Signals:')
if ~isempty(z.ulabel)
    if isstr(z.ulabel)
        disp(['   u:  ',z.ulabel])
    else
        for i=1:length(z.ulabel)
            disp(['   u',num2str(i),':  ',z.ulabel{i}])
        end
    end
end
if ~isempty(z.ylabel)
    if isstr(z.ylabel)
        disp(['   y:  ',z.ylabel])
    else
        for i=1:length(z.ylabel)
            disp(['   y',num2str(i),':  ',z.ylabel{i}])
        end
    end
end
if ~isempty(z.xlabel)
    if isstr(z.xlabel)
        disp(['   x:  ',z.xlabel])
    else
        for i=1:length(z.xlabel)
            disp(['   x',num2str(i),':  ',z.xlabel{i}])
        end
    end
end
if ~isempty(z.tlabel)
    disp(['Time:  ',z.tlabel])
end
end

function [N,ny,nu,nx]=size(z,dim)
%SIZE returns the sizes nn=[N,ny,nu,nx]
%   [N,ny,nu,nx]=size(z)
%   n=size(z,dim)  % dim=1,2,3,4

N=size(z.y,1);
ny=size(z.y,2);
nu=size(z.u,2);
nx=size(z.x,2);
if nargin==2
    if dim==1; N=N; end
    if dim==2; N=ny; end
    if dim==3; N=nu; end
    if dim==4; N=nx; end
end
if nargout==0,
   disp(['N = ',num2str(N)])
   disp(['ny = ',num2str(ny)])
   disp(['nu = ',num2str(nu)])
   disp(['nx = ',num2str(nx)])
end
end

function play(z)
%PLAY starts media player for file in field moviefile
if isempty(z.moviefile)
   error('sig.play: empty field fieldname')
end
if ~exist(z.moviefile)==2
   error(['sig.play: file ',z.moviefile,' not found'])
end
if ispc
   p=which(z.moviefile);
   if isempty(p)
     error('SIG.PLAY: animation file not found on path')
   end
   dos(['start ',p])
else
   error('SIG.PLAY: animations can only be started in Windows')
end


end

function h=plot(varargin)
%PLOT for SIG objects
%   Generates a matrix with plots of u(i) versus y(j)
%   Use yplot for just plotting y
%   See sig.sigplot for help
h=sigplot(varargin{:});
end

function h=stem(varargin)
%STEM calls sig.sigplot with view set to stem
%   See sig.sigplot for help
h=sigplot(varargin{:},'view','stem');
end

function h=staircase(varargin)
%STAIRCASE calls sig.sigplot with view set to staircase
%   See sig.sigplot for help
h=sigplot(varargin{:},'view','staircase');
end

function h=uplot(varargin)
%UPLOT calls sig.sigplot with the field y set to the signal in field u
for k=1:length(varargin)
    if isa(varargin{k},'sig')
        y=varargin{k};
        if isempty(y.u)
            error('No input in the SIG object to plot in UPLOT')
        end
        u=sig(y.u,y.t,[]);
        u.ylabel=y.ulabel;
        u.tlabel=y.tlabel;
        u.marker=y.marker;
        u.markerlabel=y.markerlabel;
        u.name=y.name;
        u.desc=y.desc;
        varargin{k}=u;
    end
end
h=sigplot(varargin{:});
end

function h=yplot(varargin)
%YPLOT plots the output y with input u removed
for k=1:length(varargin)
    if isa(varargin{k},'sig')
        y=varargin{k};
        if isnan(y.fs)
           yy=sig(y.y,y.t);
        else
           yy=sig(y.y,y.fs);
        end
        yy.ylabel=y.ylabel;
        yy.tlabel=y.tlabel;
        yy.marker=y.marker;
        yy.markerlabel=y.markerlabel;
        yy.name=y.name;
        yy.desc=y.desc;
        varargin{k}=yy;
    end
end
h=sigplot(varargin{:});
end

function h=xplot(varargin);
%XPLOT calls sig.sigplot with the field y set to the signal in field x
%   xplot(x1,x2,x3,...)
%
%   Essentially the same as plot(x2y(x1),x2y(x2),x2y(x3),...)

oldhold = ishold(gca);
m=0;
optarg={};
for k=1:length(varargin)
    if isa(varargin{k},'sig')
        y=varargin{k};
        if isempty(y.x)
            error('No states in the SIG object to plot in XPLOT')
        end
        %        x=x2y(y);
        x=sig(y.x,y.fs,[],y.x,y.Py,y.Px);
        x.ylabel = y.xlabel;
        varargin{k}=x;
    else
        m=m+1;
        optarg{m}=varargin{k};
    end
end
h=sigplot(varargin{:},optarg{:});
if oldhold==1
   hold(gca,'on')
else
   hold(gca,'off')
end
if nargout == 0
  clear h;
end
end

function h=xplot2(varargin)
%XPLOT2 plots a trajectory of x(i) versus x(j) over time
%   xplot2(x1,x2,...,Property1,Value1,...,[i j])
%   Note: index vector comes last
%
%   Property    Value/{Default}        Description
%   ----------------------------------------------------------------------
%   tlabel      {10}       '           Text label for tlabel time instants
%   view        {'plot'}|'plot3'       Plot method
%   scatter     {'off'}|'on'           Plot method
%   See plot for other property value pairs

h=[];
N=0; SIG=[]; K=0; optvar='';
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
opt=struct('axis',0,'scatter','off','conf',0,'Xlim',[],'Ylim',[],'col','bgrmcyk','linewidth',[],'fontsize',[],'view','plot','legend',[],'tlabel',10);
opt=optset(opt,optvar);
if ~isnumericscalar(opt.tlabel)
   error('sig.xplot2: tlabel must be a real scalar')
end
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

if isnumeric(varargin{end}) & sum(size(varargin{end}))==3
   ind=varargin{end};
else
   ind=[1 2];
end
i=ind(1);
j=ind(2);
col=opt.col;
xmin=Inf; xmax=-Inf; ymin=Inf; ymax=-Inf;
for k=1:length(varargin)
    kk=1+ceil(rem(k-1,length(opt.col)));
    if isa(varargin{k},'sig')
        y=varargin{k};
        if isempty(y.x)
            error('No states in the SIG object to plot in XPLOT')
        end
        if strcmp(opt.view,'plot')
	       hh=plot(y.x(:,i),y.x(:,j),[col(kk),'.-'],'parent',opt.axis,'linewidth',opt.linewidth);
           h(k)=hh(1);
        elseif strcmp(opt.view,'plot3')
	       plot3(y.x(:,i),y.x(:,j),zeros(size(y.x,1),1),[col(kk),'.-'],'parent',opt.axis,'linewidth',opt.linewidth)
	    else
            error(['sig.plot: unknown option for view: ',opt.view])
        end
        hold(opt.axis,'on')
        if opt.tlabel>0
            for t=1:round(size(y)/opt.tlabel):size(y)
               text(y.x(t,i),y.x(t,j),num2str(y.t(t)),'parent',opt.axis,'color',col(kk),'fontsize',opt.fontsize)
            end
        end
        ymax=max([ymax y.x(:,j)']);
        ymin=min([ymin y.x(:,j)']);
        xmax=max([xmax y.x(:,i)']);
        xmin=min([xmin y.x(:,i)']);
        if ~isempty(y.Px)
            if opt.conf>0 & opt.conf<=100
                confellipse2(y.x(:,[i j]),y.Px(:,[i j],[i j]),col(kk),opt.axis,opt.conf,opt.linewidth);
                a=axis;
                xmin=min([xmin a(1)]); xmax=max([xmax a(2)]);
                ymin=min([ymin a(3)]); ymax=max([ymax a(4)]);
            end
        elseif ~isempty(y.xMC)
            if strcmp(opt.scatter,'on')  % scatter always interpolated
                xii=y.xMC(:,:,i);
                xjj=y.xMC(:,:,j);
                plot(xii(:),xjj(:),['.',opt.col(kk)],'parent',opt.axis);
            end
            if opt.conf>0 & opt.conf<=100
                N=size(y.xMC,2);
                for tt=1:N
                    Px(tt,:,:)=cov(squeeze(y.xMC(:,tt,[i j])));
                end
                confellipse2(y.x(:,[i j]),Px,col(kk),opt.axis,opt.conf,opt.linewidth);
            end
        end
        if ~isempty(y.xlabel)
	    xlabel(gca,y.xlabel{i})
            ylabel(gca,y.xlabel{j})
        end
    end
end


if oldhold==1
     hold(opt.axis,'on')
else
     try
         axis(opt.axis,[xmin-0.1*(xmax-xmin) xmax+0.1*(xmax-xmin) ymin-0.1*(ymax-ymin) ymax+0.1*(ymax-ymin)])
     end
     hold(opt.axis,'off')
end
if nargout==0
  clear h;
end
end


%==================================
%-----Low-level data operations----
%==================================

function z=uplus(arg)
%UPLUS is the unitary plus
z=arg;
end

function z=uminus(arg)
%UMINUS is the unitary minus
z=arg;
z.y=-z.y;
z.yMC=-z.yMC;
end

function z=power(arg1,arg2)
%POWER
z=arg1;
z.y=power(z.y,arg2);
end

function z=mpower(arg1,arg2)
%MPOWER
z=arg1;
z.y=mpower(z.y,arg2);

end

function z=plus(arg1,arg2)
%PLUS adds two signals, a signal and a constant or a signal and a distribution
%  z1+z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1+k  if k is a numeric scalar
%  z1+e  if e is a vector with length(e)==N1
%  z1+E  if E is a PDF object, noise is added to y and yMC in which case
%        Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('plus',arg1,arg2);

end

function z=minus(arg1,arg2)
%MINUS subtracts two signals, a signal and a constant or a signal and a dist
%  z1-z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1-k  if k is a numeric scalar
%  z1-e  if e is a vector with length(e)==N1
%  z1-E  if E is a PDF object, noise is added to y and yMC in which case
%        Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('minus',arg1,arg2);

end

function z=times(arg1,arg2)
%TIMES multiplies two signals, a signal and a constant or a signal and a dist
%  z1.*z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1.*k  if k is a numeric scalar
%  z1.*e  if e is a vector with length(e)==N1
%  z1.*E  if E is a PDF object, noise is added to y and yMC in which case
%         Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('times',arg1,arg2);
end

function z=mtimes(arg1,arg2)
%MTIMES multiplies two signals, a signal and a constant or a signal and a dist
%  z1.*z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1.*k  if k is a numeric scalar
%  z1.*e  if e is a vector with length(e)==N1
%  z1.*E  if E is a PDF object, noise is added to y and yMC in which case
%         Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('mtimes',arg1,arg2);

end

function z=rdivide(arg1,arg2)
%RDIVIDE divides two signals, a signal and a constant or a signal and a dist
%  z1./z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1./k  if k is a numeric scalar
%  z1./e  if e is a vector with length(e)==N1
%  z1./E  if E is a PDF object, noise is added to y and yMC in which case
%         Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('rdivide',arg1,arg2);

end

function z=divide(arg1,arg2)
%DIVIDE is the same as RDIVIDE for SIG objects
%  z1/z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1/k  if k is a numeric scalar
%  z1/e  if e is a vector with length(e)==N1
%  z1/E  if E is a PDF object, noise is added to y and yMC in which case
%         Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('rdivide',arg1,arg2);
end

function z=mrdivide(arg1,arg2)
%MRDIVIDE is the same as MRDIVIDE for SIG objects
%  z1/z2 if z2 is a SIG object with N1==N2, fs1==fs2, t1==t2
%  z1/k  if k is a numeric scalar
%  z1/e  if e is a vector with length(e)==N1
%  z1/E  if E is a PDF object, noise is added to y and yMC in which case
%         Monte Carlo simulations are done
%  Both fields y and yMC are affected by these operations

z=evalfun2('rdivide',arg1,arg2);
end

function z=evalfun2(op,arg1,arg2)
%EVALFUN2 contains common code for plus, minus, times and divide

if isa(arg1,'sig')
    z1=arg1;
    z2=arg2;
else
    z1=arg2;
    z2=arg1;
end
if isnumericscalar(z2)  % op on a constant
    z=z1;
    z.y=feval(op,z.y,z2);
    z.yMC=feval(op,z.yMC,z2);
elseif isnumeric(z2)    % op on a vector
    if size(z1.y)==size(z2)
        z=z1;
        z.y=feval(op,z.y,z2);
        if z.MC>0 & ~isempty(z.yMC);
            tmp(1,:,:)=z2;
            z2MC=repmat(tmp,[z.MC 1 1]);
            z.yMC=feval(op,z.yMC,z2MC);
        end
    else
        error(['Using ',upper(op),' on a vector/matrix and a SIG object, size(z.y) must coincide with the size of the vector/matrix'])
   end
elseif isa(z2,'sig') % op on two signals
    if all(size(z1.y)==size(z2.y)) && all(size(z1.fs)==size(z2.fs)) && numel(z1.t)==numel(z2.t)
        z=z1;
        z.y=feval(op,z.y,z2.y);
%        z.x=feval(op,z.x,z2.x);
%        z.Px=z.Px+z2.Px;
%        z.Py=z.Py+z2.Py;
        if z1.MC==z2.MC
           z.yMC=feval(op, z1.yMC, z2.yMC);
        elseif z1.MC>0 && z2.MC>0 && ~isempty(z1.yMC) && ~isempty(z2.yMC)
           if z1.MC~=z2.MC
               MC=min([z1.MC,z2.MC]);
               z.yMC=feval(op,z1.yMC(1:MC,:,:),z2.yMC(1:MC,:,:));
           else
               z.yMC=feval(op,z1.yMC,z2.yMC);
           end
        elseif z1.MC==0 && ~isempty(z2.yMC);
            MC=z2.MC;
            tmp(1,:,:)=z1.y;
            z1MC=repmat(tmp,[MC 1 1]);
            z.yMC=feval(op,z2.yMC,z1MC);
        elseif z2.MC==0 && ~isempty(z1.yMC);
            MC=z1.MC;
            tmp(1,:,:)=z2.y;
            z2MC=repmat(tmp,[MC 1 1]);
            z.yMC=feval(op,z1.yMC,z2MC);
        else
            error('SIG: Check that the MC number coincides with the actual number of MC simulations in the yMC field')
        end
    else
        if size(z1.y,1)~=size(z2.y,1)
           error(['SIG: Using ',upper(op),' on two SIG objects, assure that the lengths are equal (N1=',num2str(size(z1.y,1)),', N2=',num2str(size(z2.y,1)),').'])
        elseif size(z1.y,2)~=size(z2.y,2)
           error(['SIG: Using ',upper(op),' on two SIG objects, assure that the signal dimensions are equal (ny1=',num2str(size(z1.y,2)),', ny2=',num2str(size(z2.y,2)),').'])
        elseif size(z1.fs)~=size(z2.fs)
           error(['SIG: Using ',upper(op),' on two SIG objects, assure that the sampling frequencies are equal (fs1=',num2str(z1.fs),', fs2=',num2str(z2.fs),').'])
        else
           error(['SIG: Using ',upper(op),' on two SIG objects, assure that the time vectors are equal.'])
        end
    end
elseif isa(z2,'pdfclass') % op on a sig and a noise distribution
   z=z1;
   [N,ny,nu,nx]=size(z1);
   ny2=length(z2);
   if ny2~=ny
      error('Signal dimension ny is not equal to the noise dimension length(z2)')
   end
   e=rand(z2,N);
   z.y=feval(op,z.y,e);
   if z.MC==0;
      z.yMC=[];
   else
      if isempty(z.yMC)
         tmp(1,:,:)=z.y;
         z.yMC=repmat(tmp,[z.MC 1 1]);
      end
      eMC=rand(z2,z.MC*N);
      eMC=reshape(eMC(:),[z.MC N ny]);
      z.yMC=feval(op,z.yMC,eMC);
   end
else
   error(['Inappropriate argument for ',upper(op),' in SIG objects'])
end
end


%==============================================
%-----OPERATORS FOR STOCHASTIC OBJECTS----------
%==============================================

function zr=rand(z,n)
%RAND returns a cell array of signal realizations taken from MC data
%   zr=rand(z,n);
%   n is the number of realizations (default n=1)
%   zr is a cell array of length n

if nargin<2; n=1; end
if isempty(z.yMC)
   error('No MC data in input argument')
end
MC=size(z.yMC,1);
if isnan(z.fs)
  in2=z.t;
else
  in2=fs;
end
for k=1:n
   ind=ceil(MC*rand(1,1));
   if isempty(z.xMC)
      x=[];
   else
      x=shiftdim(z.xMC(ind,:,:),1);
   end
   zr{k}=sig(shiftdim(z.yMC(ind,:,:),1),in2,z.u,x);
   zr{k}=inherit(zr{k},z);
end


end

function z2=fix(z1)
%FIX removes the Monte Carlo data from the SIG object
z2=z1;
z2.yMC=[];
z2.MC=0;
end

function z2=mean(z1)
%MEAN returns the mean signal of the Monte Carlo data
%z2=z1;
if isempty(z1.yMC)
   error('No MC data in input argument')
end
z2=sig(shiftdim(mean(z1.yMC,1),1),z1.fs,z1.u,z1.x,[]);
%z2.y=shiftdim(mean(z2.yMC,1),1);
%z2.yMC=[];
z2.MC=0;
end

function z2=E(z1)
%E returns the mean signal of the Monte Carlo data
z2=mean(z1);
end

function z2=std(z1)
%STD returns the standard deviation of the Monte Carlo data
if isempty(z1.yMC)
   error('No MC data in input argument')
end
z2=sig(shiftdim(std(z1.yMC,0,1),1),z1.fs,z1.u,z1.x,[]);
z2.MC=0;
end

function z2=var(z1)
%VAR returns the variance of the Monte Carlo data
if isempty(z1.yMC)
   error('No MC data in input argument')
end
z2=sig(shiftdim(std(z1.yMC,0,1),1).^2,z1.fs,z1.u,z1.x,[]);
z2.MC=0;
end


%==================================
%----High-level data operations----
%==================================

function z2=sample(z1,fs,varargin)
%SAMPLE interpolates y(t) to y[kT]
%   This function makes a call to interp
%   The help interp for property value pairs

if ~isnumericscalar(fs)
   error('fs must be a numeric scalar')
end
N=ceil((z1.t(end)-z1.t(1))*fs);
t2=z1.t(1)+(0:N-1)/fs;
z2=interp(z1,t2,varargin{:});
z2.fs=fs;  % Set this explicitely
z2=inherit(z2,z1,'Sampled data');

end

function z2=interp(z1,t2,varargin)
%INTERP interpolates y(t1) to y(t2)
%   z2=interp(z1,t2,Property1,Value1,...)
%
%   Interpolation is based on either a band-limited assumption, where
%   perfect reconstruction and re-sampling can be done, a spline
%   interpolation, or using an assumption of intersample behaviour. This
%   can be zero-order hold for piece-wise constant signal, first-order hold
%   for piece-wise linear signal.
%
%   Property   Value       Description
%   ---------------------------------------------------------
%   method    'BL' | {'hold'} | 'spline'
%   degree    {0} | 1      Degree in hold function
%
%   Examples:
%     t1=sort(rand(20,1));  % Non-uniformly sampled data
%     y1=sig(sin(2*pi*t1),t1);
%     t2=0.1:0.01:0.9;
%     y2FOH=interp(y1,t2,'method','hold','degree',1);
%     y2spline=interp(y1,t2,'method','spline');
%     plot(y1,y2FOH,y2spline)
%
%   See also: resample, decimate

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $


[N,ny,nu,nx]=size(z1);
if isempty(z1.t)
   t1=(1:N)/z1.fs;
else
   t1=z1.t;
end
if ~isnumeric(t2)
   error('t2 must be numeric')
end
if 0 & (min(t2)<min(t1)-eps | max(t2)>max(t1)+eps)  %XXX
     disp(['min(t1) = ',num2str(min(t1))])
     disp(['max(t1) = ',num2str(max(t1))])
     disp(['min(t2) = ',num2str(min(t2))])
     disp(['max(t2) = ',num2str(max(t2))])
     disp('For interpolation, the vector t2 must have all elements in the interval specified by t1')
end

opt=struct('method','hold','degree',0);
opt=optset(opt,varargin);

w=[z1.y z1.u z1.x];
nw=size(w,2);
for k=1:nw  % Loop over all signals
    y1=w(:,k);  % Vector signal

  if strcmpi(opt.method,'bl')
    if any(abs(diff(diff(t1)))>1e-12) % Non-uniform sampling
        % Use the relation y1(t1)=sum_k y(kT)sinc((t1-kT)/T)
        %                  y2(t2)=sum_k y(kT)sinc((t2-kT)/T)
        N=length(y1);
        kT=t1(1):(t1(end)-t1(1))/(N-1):t1(end);
        T=kT(2)-kT(1);
        tm=repmat(t1,1,N)-repmat(kT,N,1);
        ind0=find(tm==0);
        ind=find(tm~=0);
        tm(ind0)=tm(ind0)+eps;
        A=sin(pi*tm/T)./(pi*tm/T);
        %ind=find(isnan(A));
        A(ind0)=1;
        ykT=A\y1;
        tm=repmat(reshape(t2,[length(t2) 1]),1,N)-repmat(kT,length(t2),1);
        A=sin(pi*tm/T)./(pi*tm/T);
        ind=find(isnan(A));
        A(ind)=1;
        y2=A*ykT;
    else     % Dedicated algorithms for uniform sampling
        % Use the relation y2(t2)=sum_k y1(kT)sinc((t2-kT)/T)
        N=length(y1);
        kT=t1(:)';
        T=kT(2)-kT(1);
        tm=repmat(reshape(t2,[length(t2) 1]),1,N)-repmat(kT,length(t2),1);
        A=sin(pi*tm/T)./(pi*tm/T);
        ind=find(isnan(A));
        A(ind)=1;
        y2=A*y1;
    end
  elseif strcmpi(opt.method,'hold')
    n=opt.degree;
    if 1%any(abs(diff(diff(t1)))>1e-12) % Non-uniform sampling
        if n==0; % Zero order hold
            for i=1:length(t2);
                [dum,ind]=min(abs(t2(i)-t1));
                y2(i,:)=y1(ind,:);
            end
        elseif n==1  % First order hold
            for i=1:length(t2);
                ind=find(t2(i)<t1);
%if isempty(ind), t1,t2(i), end
                if isempty(ind)
                   y2(i,:)=y1(end,:);
                elseif ind(1)==1
                   y2(i,:)=y1(1,:);
                else
                   ind2=ind(1);
                   ind1=ind2-1;
                   tau=t1(ind2)-t1(ind1);
                   y2(i,:)=y1(ind1,:)+(y1(ind2,:)-y1(ind1,:))*(t2(i)-t1(ind1))/tau;
                end
            end
        end
    else     % Dedicated algorithms for uniform sampling
    end
  elseif strcmpi(opt.method,'spline')
    y2=interp1(t1,y1,t2,'spline','extrap');
  else
     error(['Unknown method ',method])
  end
  wip(:,k)=y2;
end
y=wip(:,1:ny);
u=wip(:,ny+1:ny+nu);
x=wip(:,ny+nu+1:ny+nu+nx);
if ~isempty(z1.yMC)
    for k=1:z1.MC
        if ~isempty(z1.xMC)
             ztmp=interp( sig(shiftdim(z1.yMC(k,:,:),1),z1.t,[],shiftdim(z1.xMC(k,:,:),1)),t2 ,varargin{:});
        else
             ztmp=interp( sig(shiftdim(z1.yMC(k,:,:),1),z1.t),t2 ,varargin{:});
        end
        yMC(k,:,:)=ztmp.y;
        xMC(k,:,:)=ztmp.x;
    end
else
    yMC=[]; xMC=[];
end
z2=sig(y,t2,u,x,yMC,xMC);
z2=inherit(z2,z1,'Interpolated data');

end

function [zd,ztrend,lsfit]=detrend(z,order);
%DETREND removes trends in non-stationary time series
%   [zd,ztrend,lsfit]=detrend(y,order);
%   A polynomial model of a certain order is estimated
%   by the least squares method and subtracted from data.
%
%   The polynomial p(t) in ztrend that is fitted to z is in absolute time,
%   not in sample index, so the function works also for
%   non-uniformly sampled data.
%   Note: detrending is only applied to the output data y, not the input or state.
%   The detrending is applied to each output signal separately.
%
%   z      Input data as a SIG object
%   order  Order of polynomial, 0 (default) for subtracting mean,
%          1 for subtracting linear trend and so on
%   zd     Detrended data in SIG object
%   ztrend The estimated trend as a SIG object
%   lsfit  Least squares loss function
%
%   Examples:
%     load genera
%     [yd2,trend2,lsfit2]=detrend(y1,2);
%     [yd3,trend3,lsfit3]=detrend(y1,3);
%     [lsfit2,lsfit3]
%     plot(y1,yd3,trend2,trend3,'Ylim',[-1000 4000])


y=z.y;
if isempty(z.t)
   t=(1:length(y))/z.fs;
else
   t=z.t;
end

if nargin<2; order=0; end

ny=size(y,2);
N=size(y,1);

for j=0:order;
     Phi(:,j+1)=t(:).^j;
end
%Phi=t(:).^(0:order);
for k=1:ny
   theta=Phi\y(:,k);
   ytrend(:,k)=Phi*theta;
   yd(:,k)=y(:,k)-ytrend(:,k);
   lsfit(k)=yd(:,k)'*yd(:,k)/N;
end
if ~isempty(z.yMC)
    for k=1:z.MC
        [zdtmp,zttmp]=detrend( sig(shiftdim(z.yMC(k,:,:),1),1),order);
        ydMC(k,:,:)=zdtmp.y;
        ytMC(k,:,:)=zttmp.y;
    end
else
    ydMC=[]; ytMC=[];
end

zd=sig(yd,t,z.u,z.x,ydMC);
zd.fs=z.fs;
zd=inherit(zd,z,'Detrended data');
ztrend=sig(ytrend,t,z.u,z.x,ytMC);
ztrend.fs=z.fs;
ztrend=inherit(ztrend,z,'Trend in data');

end

function zw=window(z,varargin)
%WINDOW computes and applies a data window to SIG object
%   yw=window(z,type,n)
%   Windowing of data returns yw=y.*w, where w is obtained from getwindow
%   For IO systems, also xw is windowed
%
%   type is 'box', 'hanning', 'hamming','kaiser',
%           'blackman', 'bartlett' or 'spline'
%   spline uses a uniform window convolved with itself n times.
%   See getwindow for detailed help on type and n
%
%   Examples:
%   See also: getwindow

if nargin<2
  type='hanning';
end

y=z.y;
[N,ny,nu,nx]=size(z);
w=getwindow(N,varargin{:});

uw=[]; xw=[];
yw=y.*(w*ones(1,ny));
if ~isempty(z.u)
    uw=z.u.*(w*ones(1,nu));
end
if ~isempty(z.x)
    xw=z.x.*(w*ones(1,nx));
end
%xr=wr(:,ny+nu+1:ny+nu+nx);

if ~isempty(z.yMC)
    for k=1:z.MC
        if ~isempty(z.xMC)
             ztmp=window( sig(shiftdim(z.yMC(k,:,:),1),1,[],shiftdim(z.xMC(k,:,:),1)),varargin{:});
        else
             ztmp=window( sig(shiftdim(z.yMC(k,:,:),1),1),varargin{:});
        end
        ywMC(k,:,:)=ztmp.y;
        xwMC(k,:,:)=ztmp.x;
    end
else
    ywMC=[]; xwMC=[];
end
zw=sig(yw,z.t,uw,xw,ywMC,xwMC);
zw.fs=z.fs;
zw=inherit(zw,z,['Windowed']);
end

function z2=decimate(z1,n,varargin)
%DECIMATE decimates the signal a factor n, using anti-alias filtering
%   z2=decimate(z1,n)
%   This is equivalent to z2=resample(z1,n)
%   downsample also decimates a signal, but keeps all samples in
%   concatenated vectors
%   See help resample

z2=resample(z1,n);
z2=inherit(z2,z1,'Decimated data');
end

function z2=resample(z1,n,m)
%RESAMPLE resamples uniformly sampled signal using a band-limitation assumption
%   z2=resample(z1,n,m)
%
%   Resampling definition: y[k]=y(kT) to y[l]=y(l*n/m*T)
%   resample(z,n,1) decimates a factor n
%   resample(z,1,m) upsamples a factor m
%   Anti-alias filtering applied if n/m>1 automatically
%
%   Examples:
%     t1=(0:1:100)';
%     y1=sin(0.1*t1);
%     z1=sig(y1,t1);
%     z2=resample(z1,7,3);
%     plot(z1,z2)
%   See also: interp

if nargin<3, m=1; end
if ~isnumericscalar(n) | ~isnumericscalar(m) | round(n)~=n | round(m)~=m | m<1 | n<1
   error('n and m must be positive integers')
end
y=z1.y;
fs=z1.fs;
if isnan(fs)
   error('RESAMPLE applies only to discrete time signals with Uniform sampling')
end
T=1/fs;
if isempty(z1.t)
   t=(1:length(y))/fs;
else
   t=z1.t;
end
[N,ny,nu,nx]=size(z1);
w=[z1.y z1.u z1.x];
nw=size(w,2);
for k=1:nw
    u=w(:,k);
    % Upsampling
    if m>1
        U=fft(u);
        if N/2==round(N/2) % even N
           Um=[U(1:N/2+1);zeros((m-1)*N,1);U(N/2+2:N)];
        else
           Um=[U(1:(N+1)/2);zeros((m-1)*N,1);U((N+1)/2+1:N)];
        end
        u=real(ifft(Um));
        Nup=length(u);
        t=t(1): (t(end)-t(1))/(Nup-1): t(end);
    end
    % Downsampling
    if n>1
        Nup=length(u);
        %Zero pad to make length(u)/n an even integer
        u=[u;zeros(ceil(Nup/2/n)*n*2-Nup,1)];
        Nnew=length(u);
        Nn=Nnew/n;
        U=fft(u);
        Un=[U(1:Nn/2+1);U(end-Nn/2+1:end)];
        u=real(ifft(Un));
        u=u(1:Nn);
        t=t(1): (t(end)-t(1))/(Nn-1): t(end);
%        t=(0:1/m:Nnew/m)'*Nup/Nnew*T;
%        t=t(1:end-1);
%        t=t(1:n:end);
        u=u*Nnew/Nup;
    end
    wr(:,k)=(m/n)*u;
end
yr=wr(:,1:ny);
ur=wr(:,ny+1:ny+nu);
xr=wr(:,ny+nu+1:ny+nu+nx);
if ~isempty(z1.yMC)
    for k=1:z1.MC
        if ~isempty(z1.xMC)
             ztmp=resample( sig(shiftdim(z1.yMC(k,:,:),1),1,[],shiftdim(z1.xMC(k,:,:),1)),n,m);
        else
             ztmp=resample( sig(shiftdim(z1.yMC(k,:,:),1),1),n,m);
        end
        yrMC(k,:,:)=ztmp.y;
        xrMC(k,:,:)=ztmp.x;
    end
else
    yrMC=[]; xrMC=[];
end
z2=sig(yr,t,ur,xr,yrMC,xrMC);
z2.fs=fs*m/n;
z2=inherit(z2,z1,['Resampled a factor ',num2str(m),'/',num2str(n)]);
end

function z2=downsample(z1,n,varargin)
%DOWNSAMPLE downsamples the signal a factor n, keeping all samples
%   z2=downsample(z1,n)
%   The new SIG object has n*ny outputs and n*nu inputs.
%   The sampling frequency is reduced to fs/n
%   The length of the new signal is floor(N/n)
%   Compare to deciate, which throws away samples
%   See also resample, decimate

[N,ny]=size(z1.y);
ym=z1.y';
ym=ym(:);
Nm=floor(N/n);
y=reshape(ym(1:Nm*n*ny),n*ny,Nm)';
[N,nu]=size(z1.u);
if nu>0
   um=z1.u';
   um=um(:);
   Nm=floor(N/n);
   u=reshape(um(1:Nm*n*nu),n*nu,Nm)';
else
   u=zeros(Nm,0);
end
x=z1.x(n:n:n*Nm,:);
z2=sig(y,z1.fs/n,u,x);
z2=inherit(z2,z1,'Downsampled data');
end



%=========================
%----Conversions----------
%=========================

function y=u2y(u);
%U2Y extracts the input signal as a SIG object
%   y=u2y(u);
%   The SIG object (y,x,u) is mapped to (u,[],[])
%
%   Example:
%     G=rand(lss([3 1 0 2],1)); % System
%     u=getsignal('prbs');     % Input
%     z=simulate(G,10)         % Simulation -> z=(y,x,u)
%     u=u2y(z)                 % Recovered input


y=sig(u.u,u.t,[]);
y.ylabel=u.ulabel;
y.tlabel=u.tlabel;
y.marker=u.marker;
y.markerlabel=u.markerlabel;
y.name=u.name;
y.desc=u.desc;
y.fs=u.fs;
end

function y=x2y(u);
%X2Y extracts the state as a SIG object
%   y=x2y(x);
%   The SIG object (y,x,u) is mapped to (x,x,u)
%
%   Example:
%     G=rand(lss([3 1 0 2],1)); % System
%     u=getsignal('prbs');     % Input
%     z=simulate(G,10)         % Simulation -> z=(y,x,u)
%     u=u2y(z)                 % Recovered input

y=sig(u.x,u.t,u.u,u.x,u.xMC,u.xMC);
y.Py=u.Px;
y.ylabel=u.xlabel;
y.ulabel=u.ulabel;
y.tlabel=u.tlabel;
y.marker=u.marker;
y.markerlabel=u.markerlabel;
y.name=u.name;
y.desc=u.desc;
y.fs=u.fs;

end

function Y=sig2ft(z,varargin);
%SIG2FT computes Fourier Transform (approximation) of a signal
%   Y=sig2ft(y,Property1,Value1,...);
%   This function is used in the FT constructor
%
%   For discrete time signals, the function basically performs a
%   fft for each signal dimension of a zero-padded signal.
%   For continuous time signals (fs=NaN), the Dirichlet approximation
%   of the signal is computed, which basically assumes that the signal
%   consists of impulses z.y at the times specified in z.t
%
%   Y is an object with fields Y for DTFT and f for frequency
%   f is computed uniformly between 0 and fs/2
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   Nf         {8*N}      Number of grid points in f
%                         Use Nf=N to obtain the usual DFT grid
%   f                     Frequency grid (overrides the default one)
%                         This can be used to zoom in interesting parts
%
%   Examples:
%      y=sin(0:0.4:20);
%      subplot(2,1,1),  plot(abs(fft(y)))
%      subplot(2,1,2),  plot(ft(sig(y)))
%
%   See also: ft, ft.plot

%$ Revision: 28-Oct-2019 $


y=z.y;
[Ny,ny]=size(y);
fs=z.fs;
t=z.t;
if z.MC>0
    yMC=z.yMC;
else
    yMC=[];
end
Ts=1/fs;

opt=struct('Nf',8*Ny,'f',[]);
opt=optset(opt,varargin);

Nf=opt.Nf;
if Nf/2~=round(Nf/2)
   Nf=Nf+1;
   disp(['Warning: Nf must be even, and is changed to ',num2str(Nf)])
end
if ~isnan(fs)  % DTFT
    if isempty(opt.f)
        for k=1:ny
            Ytmp=fft([y(:,k);zeros(Nf-Ny,1)]);
            Y(:,k)=Ytmp(1:Nf/2+1);
            f= (0:Nf/2)'/Nf/Ts;
        end
    else
        f=opt.f(:);
        W=exp(-i*2*pi*f*t');
        Y=W*y;  % Works also when ny>1
    end
else % % Dirichlet approximation of FT
    if isempty(opt.f)
        f= (0:Nf/2)'/t(end);
    else
        f=opt.f;
    end
    W=exp(-i*2*pi*f*t(:)');
    Y=W*y;     % Works also when ny>1
end

YMC=[];
if ~isempty(z.yMC)
    for k=1:size(z.yMC,1)
        Ytmp=sig2ft(sig(shiftdim(z.yMC(k,:,:),1),fs));
        YMC(k,:,:)=Ytmp.Y;
     end
end
Y=ft(Y,f,YMC);
Y=inherit(Y,z,'SIG -> FT');

end

function c=sig2covfun(z,varargin)
%SIG2COVFUN estimates the (cross-)covariance function of (the columns in) y
%   c=sig2covfun(y,Property1,Value1,...)
%
%   z  signal object
%   c  covariance object
%   P=sig2covfun(z,'taumax',0) returns the (spatial) covariance matrix
%
%   Property   Value        Description
%   ---------------------------------------------------------
%   taumax     {30}         Maximum lag for which the covariance function is computed
%   fs         {y.fs}       Sampling frequency (overrides fs specified in SIG y)
%   MC         {100}        Number of Monte Carlo simulations to compute confidence bound
%   method     {'direct'}   Direct summation in the time domain
%              'conv'       Summation in the time domain using conv
%              'freq'       Convolution computed in the frequency domain
%
%   Examples:
%
%   See also:

fs=z.fs;
y=z.y;

opt=struct('taumax',30,'MC',100,'method','direct','fs',fs);
opt=optset(opt,varargin);

[N,ny]=size(y);
%if N<ny, disp('Warning: N<ny'), end
taumax=opt.taumax;
Rtot=zeros(taumax+1,ny,ny);
for j=1:ny
    for i=j:ny
        if strncmpi(opt.method,'direct',2);  % Time domain, direct method
            for k=0:taumax
                R(k+1)=y(1:end-k,i)'*y(1+k:end,j)/(N-k);
            end
            R=R(:);
        elseif strncmpi(opt.method,'conv',2); % Time domain, conv
            R=conv(y(:,i),y(end:-1:1,j));
            R=R(N:N+taumax)./(N-(0:taumax)');
        elseif strncmpi(opt.method,'freq',2) % frequency domain
            N2 = 2^(ceil(log2(N)));
            Yi=fft(y(:,i),N2);
            Yj=fft(y(:,j),N2);
            Phi=(Yi.*conj(Yj));
            R=real(ifft(Phi));
            R=R(1:taumax+1)./(N-(0:taumax)'); % Unbiased estimate
            %R=R(1:taumax+1)/N;
        else
            error('method is one of ''direct'', ''conv'' or ''freq''')
        end
        Rtot(:,i,j)=R;
        Rtot(:,j,i)=R;
    end
end
R=Rtot;
RMC=[];
if ~isempty(z.yMC)
    for k=1:z.MC
        Rtmp=sig2covfun(sig(shiftdim(z.yMC(k,:,:),1),fs),varargin{:});
        RMC(k,:,:)=Rtmp.R;
     end
end
c=covfun(R,0:taumax,RMC);
c=inherit(c,z,' SIG -> COVFUN');

end

function Phi=sig2spec(z,varargin)
%SIG2SPEC performs spectral analysis of the signal object z
%   Phi=sig2spec(z,Property1,Value1,...)
%
%   For both uniformly and non-uniformly sample data, the periodogram
%   and its smoothed version can be computed.
%   For uniformly sampled data, also the Welch and
%   Blackman-Tukey methods are available.
%
%   Phi is a SPEC spectrum object
%   y   is a SIG signal object
%
%   Property  Value/{Default}      Description
%   ---------------------------------------------------------------------------
%   MC        {100}                Number of Monte Carlo simulations
%   M         {min([N/5 max([N/10 30])])}
%                                  Smoothing parameter, M=1 recovers the
%                                  periodogram, larger M gives less
%           detail but better averaging. Default is length(y)/30
%           For method 5, M is the threshold mu
% method    1 or 'periodogram'  squared magnitude of Fourier transform
%           2 or 'Blackman-Tukey' smoothed version of periodogram
%             with window win (default 'hamming') with width M
%           3 or 'Welch' averaged periodogram over different signal
%             segments of length M, windowed by win (default 'hamming')
%           4 or 'smoothing' applies a smoothing window directly on the
%           periodogram
%           5 or 'cepstrum' applies thresholding to the cepstrum coefficients
% overlap   Overlap used in Welch method (default 0)
% fs        Sampling frequency, scales the frequency axle f (default fs=2).
%           Overrides the fs specified in struct y.
% win       Window used in method 2 and 3. See help window for options
%           ('hamming' default)
%
%   Uncertainty is computed in Welch method interpreting the periodogram from the
%   different segments as Monte Carlo data. For smoothed periodogram, the
%   variance over each position of the sliding window is used as uncertainty.
%
%   Example:
%
%   See also:

%$ Revision: 28-Oct-2019 $


fs=z.fs;
N=size(z,1);
if isnan(fs)
    defmethod=1;
else
    defmethod=3;
end

M=ceil(min([N/5 max([N/10 30])]));
opt=struct('MC',100,'fs',fs,'method',defmethod,'M',M,'win','hamming','overlap',0);
opt=optset(opt,varargin);

if isnan(fs)
   if opt.method==2 | opt.method==3
      error('sig.sig2spec: Welch and Blackman-Tukey are not available for non-uniform data')
   else
      fs=1/z.t(end);
   end
end

method=opt.method;
if ischar(method)
    method=lower(method);
    if findstr(method,'pe'); m=1;
    elseif findstr(method,'bl'); m=2;
    elseif findstr(method,'we'); m=3;
    elseif findstr(method,'sm'); m=4;
    elseif findstr(method,'ce'); m=5;
    else error([method,' is not a valid method'])
    end
elseif method>0 & method<6
    m=method;
else
    error([method,' is not a valid method'])
end
if m==2
    Rcov=sig2covfun(z,'taumax',opt.M);
end
ny=size(z,2);

for ky=1:ny;
for ly=1:ny;
y1=z.y(:,ky);
y2=z.y(:,ly);
if m==1  | m==4  | m==5;
    N2 = 2^(ceil(log2(N)));
    Y1=fft(y1,N2);
    Y2=fft(y2,N2);
    Y1=Y2(1:N2/2+1);
    Y2=Y2(1:N2/2+1);
    f=(0:N2/2)'/N2*fs;
    Phihathat=fs/N*real(Y1.*conj(Y2));
end
if m==1
    method='periodogram';
    Phihat(:,ky,ly)=Phihathat;
    opt.M=1;
    sigma=zeros(size(Phihat));
    PhiMC=[];
elseif m==4;
    method='smoothing';
    w=getwindow(round(N/opt.M),opt.win);
    w=w/sum(w);
    Phihat(:,ky,ly)=filtfilt(w,1,Phihathat);
    Phihat2=filtfilt(w,1,Phihathat.^2);
    sigma=sqrt(abs(Phihat2-Phihat.^2));
    s=randn(1,opt.MC);
    PhiMC(:,:,ky,ly)=(sigma*s+Phihat*ones(1,opt.MC))';
elseif m==5;
    method='cepstral thresholding';
    Phi1=Phihathat(:,ky,ly);
    logphi=log(Phi1);
    ck=ifft(logphi);
    ck(1)=ck(1)+0.577216;
    mu=opt.M;
    ind=find(abs(ck)>mu);
length(ind)
    ck(ind)=zeros(size(ind));
    Phi2=exp(fft(ck));
    alpha=Phi2\Phi1;
    Phihat(:,ky,ly)=alpha*Phi2;
    PhiMC=[];
elseif m==2;
    method='Blackman-Tukey';
    w=getwindow(2*opt.M+1,opt.win);
    R=Rcov.R(:,ky,ly);
    Rw=[R(end:-1:2); R].*w;
    f=(0:opt.M)/(2*opt.M+1)*fs;
    Phihat=abs(fft(Rw))*4;
    Phihat=Phihat(1:opt.M+1);
    PhiMC=[];
elseif m==3;
    method='Welsh';
    M=ceil(opt.M/2)*2; % Even segment length
    R=ceil(N/(M-opt.overlap)); % Number of periodograms
    w=getwindow(M,opt.win);
    %sum(w),sum(w.^2),S
    w=w/sum(w)*M;
    y1=[y1(:);zeros(M*R-N,1)]; % zero pad to length M*R
    y2=[y2(:);zeros(M*R-N,1)]; % zero pad to length M*R
    f=(0:M/2)'/M*fs;
    for i=0:R-1
        Y1=fft(w.*y1(i*(M-opt.overlap)+1:i*(M-opt.overlap)+M));
        Y2=fft(w.*y2(i*(M-opt.overlap)+1:i*(M-opt.overlap)+M));
        Y1=Y1(1:M/2+1);
        Y2=Y2(1:M/2+1);
        Phihathat(:,i+1)=fs/M*(Y1.*conj(Y2));
    end
    Phihat(:,ky,ly)=mean(Phihathat,2);
    PhiMC(:,:,ky,ly)=Phihathat';
end
end
end
% Replace Monte Carlo data from method by ensemble data
if ~isempty(z.yMC)
    for k=1:z.MC
        Phitmp=sig2spec(sig(shiftdim(z.yMC(k,:,:),1),fs),varargin{:});
        PhiMC(k,:,:,:)=Phitmp.Phi;
    end
end
Phi=spec(Phihat,f,PhiMC);
Phi=inherit(Phi,z,' SIG -> SPEC');
%Phi.method=[method];

end

function Yt=sig2tfd(z,varargin)

%SIG2TFD performs a Time-Frequency Description (TFD) of the signal y
%   tfd=sig2tfd(z,Property1,Value1,...)
%   The function computes the periodogram over segments of the signal.
%   A window is applied on each segment, and the segments may overlap.
%   The output is a struct tfd with fields E (energy), t (time) and f (frequency).
%
%   The computations are very similar to Welch spectral estimate, in fact
%   taking the average over time of tfd yields the Welch spectral estimate
%   mean(tfd.E').
%   In the same way, the smoothed signal energy is computed by mean(tfd.E).
%
%   Output parameters
%   tfd    Time Frequency Description
%   tfd.E  Energy in each time frequency bin
%   tfd.f  Frequency
%   tfd.t  Time
%
%   Required input parameters
%   z   signal object with data field y
%
%   Optional parameters
%   Property  Value/{Default}  Description
%   ------------------------------------------------------------------------------
%   S         {max(N/25,128)}  Segment length in samples. The larger S, the better
%                              frequency resolution but the worse time resolution
%   overlap   {75[%]}          Overlap of each segment in percent
%   fs        {2}              Sampling frequency, scales the frequency axis f
%   win       {'hamming'}      Data window on each segment, see window for options
%
%   Examples:
%   load bach
%   sig2tfd(y,'S',2000)
%
%   See also:
%$ Revision: 28-Oct-2019 $

y=z.y;
fs=z.fs;
if isnan(fs);
    error('TFD can only be computed for discrete time signals')
end
N=length(y);
Sdefault=max([2*round(N/50) 128]);
opt=struct('S',Sdefault,'fs',fs,'win','kaiser','overlap',90);
opt=optset(opt,varargin);

S=opt.S;
fs=opt.fs;
win=opt.win;
overlap=round(S*opt.overlap/100);

M=ceil(N/S); % Number of segments
M2=floor(N/(S-overlap)); % Number of time points
t=(1:M2)/fs*(S-overlap)+ (S-overlap)/2/fs;
y=[y(:);zeros(S*M2-N,1)]; % zero pad to length S*M2
f=(0:S/2)'/S*fs;
w=getwindow(S,win);
for i=0:M2-1
    Y=fft(w.*y(i*(S-overlap)+1:i*(S-overlap)+S));
    Y=Y(1:S/2+1);
    E(:,i+1)=fs/S*real(Y.*conj(Y));
end

Yt=tfd(E,t,f);
Yt=inherit(Yt,z,'SIG -> TFD');
end

function Y=sig2ndist(z)
%SIG2NDIST converts y(t) and Pyy(t) to Y=ndist(y(t),Pyy(t))

end

function h=sigplot(varargin)
%SIGPLOT plots a signal struct
%   sigplot(z1,z2,...,Property1,Value1,...)
%
%   Property   Value       Description
%   ---------------------------------------------------------
%   view      'staircase' | {'interp'} | 'stem'
%                          Type of plot for sampled signals
%                          For continuous time signals, view is 'interp'
%   interval  {1:N}        Time interval for focus
%   axis      {gca}        Axis handle where plot is added
%   conf      [{0},100]    Confidence level from MC, 0 means no levels plotted
%   conftype  1,{2}        Confidence interval: 1 band, 2 lines
%   scatter   'on'|{'off'} Scatter plot of MC data
%   col       {'bgrmyk'}   Colors in order of appearance in z1,z2,...
%   ind       {':'}        Indeces of the signal
%   Xlim                   Limits on x axis
%   Ylim                   Limits on y axis
%   linewidth {}           Line width (default in global SIGNAL)
%   fontsize  {}           Font size  (default in global SIGNAL)
%   legend    {''}         Legend text. Default is the name fields of the signals.
%
%   Examples:
%     s1=getsignal('sin1',100);
%     subplot(2,1,1),  sigplot(s1,'view','stem')
%     s2=getsignal('sin2',100);
%     subplot(2,1,2),  sigplot(s2,'view','staircase')
%
%   See also: xplot
%

%$ Revision: 28-Oct-2019 $

N=0; SIG=[]; K=0; optvar='';
k=0;
h=[];
while k<length(varargin)
    k=k+1;
    if isstr(varargin{k})  % Property Value Pair
        K=K+1;
        optvar{K}=varargin{k};
        K=K+1;
        k=k+1;
        try
           optvar{K}=varargin{k};
        catch
	   error('Non-sig object arguments must appear in property-value pairs')
        end
    else   % vector or struct
        N=N+1;
        Z{N}=varargin{k};
    end
end
opt=struct('axis',0,'scatter','off','conf',0,'conftype',2,'Xlim',[],'Ylim',[],'col','bgrmyk','linewidth',[],'fontsize',[],'view','interp','legend',[],'ind',':');
opt=optset(opt,optvar);
if opt.axis==0; opt.axis=gca; end
oldhold = ishold(opt.axis);

Ncol=length(opt.col);
if N>Ncol  % repeat the color pattern
   opt.col=char(kron(ones(1,ceil(N/Ncol)),opt.col));
end

legon=0;
leg='';
for k=1:N;
  kk=1+ceil(rem(k-1,length(opt.col)));
  z=Z{k};
  if ~isa(z,'sig')
     error('SIG: input to plot functions must all be SIG objects')
  end
    z=z(:,opt.ind); %xxx

  [Ny,ny(k),nu(k)]=size(z);
  if k>1 % check dimensions
    if ny(k)~=ny(k-1)
      error('Signals have different output dimension')
    end
    if nu(k)~=nu(k-1)
      %error('Signals have different input dimension')
      % Not important
    end
  end
  if ny(k)>8
    error('sig.plot can only plot up to eight subplot rows, try to reduce signal dimension for y with plot(y(:,1:8)) or similar')
  end
  if nu(k)>8
    error('sig.plot can only plot up to eight subplot columns, try to reduce signal dimension for u with plot(y(:,:,1:8)) or similar')
  end

  if ~isempty(z.tlabel)
      tlabel=z.tlabel;
  else
      tlabel='Time';
  end
  for j=1:ny(k)
    nuk=max([nu(k) 1]);
    for i=1:nuk
      if nuk>1 | ny(k)>1
        opt.axis=subplot(ny(k),nuk,(j-1)*nuk+i);
        oldhold = ishold(opt.axis);
      end
      t=z.t; y=z.y(:,j);
      if ~isempty(z.yMC)
          yMC=z.yMC(:,:,j);
      end
      if nu(k)>0
        u=z.u(:,i);
      else % No input
        u=0;
      end
      if z.fs>0
        if strncmpi(opt.view,'staircase',5)
            t=kron(t(:),[1;1]); t(1)=[];
            u=kron(u(:),[1;1]); u(end)=[];
            y=kron(y(:),[1;1]); y(end)=[]; %or y(end)=[]; !
        elseif strncmpi(opt.view,'stem',4)
            t=kron(t(:)',[1;1]);
            u=kron(u(:)',[0;1]);
            y=kron(y(:)',[0;1]);
        else % interp
            % Nothing needs to be done
        end
      else  % Continuous time signal, look for impulses
        ind1=find(diff(t)==0);
        ind2=find(diff(ind1)==1);
        indimp=ind1(ind2+1);
      end
      % Start plotting
      p1=plot(t,y,['-',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth);
      h(k)=p1(1);
      set(opt.axis,'fontsize',opt.fontsize);
      xlabel(opt.axis, tlabel);
      if ~isempty(z.ylabel)
        try
          ylabel(opt.axis, z.ylabel{j});
        catch
          ylabel(opt.axis, z.ylabel);
        end
      end
      hold(opt.axis,'on')

      if isnan(z.fs) & length(u)<2  % Cannot have impulses in u and y both
         plot(t(indimp),y(indimp),['s',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth);
      end
      if strncmpi(opt.view,'stem',4)
        p1o=plot(t(2,:),y(2,:),['o',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth);
      end
      if length(u)>3
        p2=plot(t,u,['--',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth);
        if isnan(z.fs)
           plot(t(indimp),u(indimp),['s',opt.col(kk)],'parent',opt.axis,'linewidth',opt.linewidth);
        end
      end
      % MC
      if ~isempty(z.yMC)
        if strcmp(opt.scatter,'on')  % scatter always interpolated
            plot(z.t,z.yMC(:,:,j),['--',opt.col(kk)],'parent',opt.axis,'linewidth',0.5*opt.linewidth);
        end
        if opt.conf>0 & opt.conf<=100
            p=confband(z.t,z.yMC(:,:,j)','plot',opt.col(kk),opt.axis,opt.conf,opt.conftype,opt.linewidth);
            %if legon==0; legend([p1(2) p2 p4],'Nominal','Confidence bound','Median','parent',opt.axis); legon=1; end
        end
      end
      if ~isempty(z.Py)
        if opt.conf>0 & opt.conf<=100
            p=confband2(z.t,z.y(:,j),z.Py(:,j,j),'plot',opt.col(kk),opt.axis,opt.conf,opt.conftype,opt.linewidth);
        end
      end
      if nu(k)==0
        try
           title(opt.axis, z.name{j})
        end
      elseif nu(k)>1 | ny(k)>1
        try
           title(opt.axis,[z.ulabel{i},' to ',z.ylabel{j}])
        end
      end
      if k==N  % last signal
         drawnow
         ylim=get(opt.axis,'Ylim');
%         set(opt.axis,'Ylim',[ylim(1)-0.1*diff(ylim) ylim(2)+0.1*diff(ylim)])
         if ~isempty(opt.Xlim),
            set(opt.axis,'Xlim',opt.Xlim)
         end
         if ~isempty(opt.Ylim),
            set(opt.axis,'Ylim',opt.Ylim)
         end
         % markers
         if  ~isempty(z.marker)
           drawnow
           a=axis;
           plot([1;1]/z.fs*z.marker(:)',[a(3:4)]'*ones(1,length(z.marker)),'-o','linewidth',2*opt.linewidth,'parent',opt.axis)
         end

         if oldhold==1
           hold(opt.axis,'on')
         else
           hold(opt.axis,'off')
         end

      end  % end last signal case
    end  % end u loop
  end % end y loop
end  % end k loop

if ~isempty(opt.legend),
    legend(opt.axis,opt.legend{:})
elseif ~isempty(leg),
    %legend(opt.axis,leg{:})
end
if ~isempty(z.name) & ny==1 & nu<2
      title(opt.axis,z.name)
end
end


end %methods
end %class
