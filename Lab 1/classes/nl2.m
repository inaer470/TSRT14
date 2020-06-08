classdef nl2
% nl2 is a development version for implicit noise
%
%NL2 is the model object for Non-Linear time-varying systems
%
%   help nl2.nl2    gives help for the constructor
%   methods nl2    lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
   f=[];h;nn;pv;pe;px0;x0;fs;th;P;I;H=[];
   xlabel;thlabel;ulabel;ylabel;tlabel;name;desc;
end

methods

function m=nl2(f,h,nn,fs)
%NL2 constructor non-linear dynamic models
%   m=nl2(f,h,nn,fs)
%   The model is defined as
%       x+(t) = f(t,x(t),v(t),u(t);th),  v(t)~pv
%        y(t) = h(t,x(t),e(t),u(t);th),  e(t)~pe
%     E[x(0)] = x0, x(0)~px0
%   where x+(t)=x'(t) for continuous time and x+(t)=x(t+1) for discrete
%   time systems, respectively.
%   f and h are inline objects or strings with arguments
%   t,x,v,u,th, and t,x,e,u,th, respectively,
%   and nn=[nx nv nu ny nth] gives the dimensions of x,v,u,y and the parameters th.
%   Convention of dimensions:
%     x (nx,1) vector or (nx,N) matrix
%     v (nv,1) vector or (nv,N) matrix
%     u (nu,1) vector or (nu,N) matrix
%     y (ny,1) vector or (ny,N) matrix
%     e (ny,1) vector or (ny,N) matrix
%     th vector of length nth
%   Special cases:
%     nx=0, gives a signal model without state dynamics
%     ny=0, gives a dynamic system without explicit measurements
%
%   Distributions pv, pe, px0 are specified as pdfclass objects, or
%   as a covariance matrix of suitable dimension, in which
%   case a ndist object is created.
%   Other fields: name, xlabel thlabel ulabel ylabel
%
%   Examples:
%     m=nl2('-th(1)*x(1,:).^2-th(2)*x(1,:)+v(1,:)','x',[1 1 0 1 2])
%   is equivalent to
%     f=inline('-th(1)*x^2-th(2)*x+v','t','x','v','u','th');
%     h=inline('x+e','t','x','e','u','th');
%     m=nl2(f,h,[1 1 0 1 2])
%   Complete with
%     m.fs=1; m.th=[0.5 0.1]; m.x0=1; m.pv=1; m.pe=1;


%#
%#


if nargin<1
   return  % empty model
elseif isa(f,'nl2')
   m.f=f.f; m.h=f.h; m.nn=f.nn; m.fs=f.fs;
   m.H=f.H;
   m.pv=f.pv;
   m.pe=f.pe;
   m.x0=f.x0;
   m.th=f.th;
   m.px0=f.px0;
   m.xlabel=f.xlabel;
   m.ylabel=f.ylabel;
   m.ulabel=f.ulabel;
   m.name=f.name;
   m.desc=f.desc;
   return
elseif isa(f,'lss')
   mss=lss2nl2(f);
   m.f=mss.f; m.h=mss.h; m.nn=mss.nn; m.fs=mss.fs;
   m.H=mss.H;
   m.pv=mss.pv;
   m.pe=mss.pe;
   m.x0=mss.x0;
   m.px0=mss.px0;
   m.xlabel=mss.xlabel;
   m.ylabel=mss.ylabel;
   m.ulabel=mss.ulabel;
   m.name=mss.name;
   m.desc=mss.desc;
   return
elseif isa(f,'ltf')
   mss=lss2nl2(lss(f));
   m.f=mss.f; m.h=mss.h; m.nn=mss.nn; m.fs=mss.fs;
   m.H=mss.H;
   m.pv=mss.pv;
   m.pe=mss.pe;
   m.x0=mss.x0;
   m.px0=mss.px0;
   m.xlabel=mss.xlabel;
   m.ylabel=mss.ylabel;
   m.ulabel=mss.ulabel;
   m.name=mss.name;
   m.desc=mss.desc;
   return
end
if nargin<3
   error('NL2 constructor: three input arguments required')
end
if nargin<4
   fs=NaN;
end

if ~isvector(nn) | length(nn)~=5
   error('NL2 constructor: nn must be a vector of length 5')
end
nx=nn(1); nv=nn(2); nu=nn(3); ny=nn(4); nth=nn(5);
if isa(f,'inline')
   % do nothing
elseif isa(f, 'function_handle')
   % do nothing
elseif isstr(f)
   if exist(f,'file')  % m-file
      f=[f,'(t,x,v,u,th)'];
   else
       try
          f=inline(f,'t','x','v','u','th');
       catch
          error('NL2 constructor: a string argument for f must be a valid m-file or a valid inline input argument')
       end
   end
% $$$    try
% $$$       f=inline(f,'t','x','v','u','th');
% $$$    catch
% $$$        f
% $$$       error('NL2 constructor: a string argument for f must be a valid m-file or a valid inline input argument')
% $$$    end
elseif isempty(f) & nx==0
   % Static model is OK
else
   error('NL2 constructor: f must be an inline object, a string, or an m-file')
end

if isa(h,'inline')
   % do nothing
elseif isa(h, 'function_handle')
   % do nothing
elseif isstr(h)
   if exist(h,'file')  % m-file
       h=[h,'(t,x,e,u,th)'];
   else
      try
         h=inline(h,'t','x','e','u','th');
      catch
         error('NL2 constructor: a string argument for h must be a valid m-file or a valid inline input argument')
      end
   end
elseif isempty(h)
   h=[];
else
   error('NL2 constructor: h must be an inline object, a string, or an m-file')
end
x0=zeros(nx,1);
th=zeros(nth,1);
v=zeros(nv,1);
x=zeros(nx,1);
if nx>0
    try
       fdum=feval(f,0,1e-5*randn(size(x0)),1e-5*randn(nv,1),1e-5*randn(nu,1),1e-5*randn(size(th)));
%       fdum=feval(f,0,x0,zeros(nu,1),th);
    catch ME
       disp(['NL2 error code: ',ME.message])
       error('NL2 constructor: sizes in nn not consistent with f')
    end
   if size(fdum,1)~=nx
      disp(['size(f): ',num2str(size(fdum))])
      disp(['size(x,1): ',num2str(nx)])
      error('NL2 constructor: size of f not consistent with nn')
   end
   if size(fdum,2)~=1
      error('NL2 constructor: number of columns in f is not one')
   end
   try
       fdum=feval(f,[0 1],[x0 x0],zeros(nv,2),zeros(nu,2),th);
       if size(fdum,2)~=2
           error('NL2 constructor: size of f not consistent with nn')
       end
    catch ME
       disp('NL2 warning: Try to vectorize f for increased speed')
       disp(['Hint: ',ME.message])
   end
end

if ny>0
   try
      hdum=feval(h,0,1e-5*randn(nx,1),1e-5*randn(ny,1),1e-5*randn(nu,1),1e-5*randn(nth,1));
    catch ME
       disp(['NL2 error code: ',ME.message])
      error('NL2 constructor: sizes in nn not consistent with h or error in evaluating h')
   end
   if size(hdum,1)~=ny
      disp(['size(h): ',num2str(size(hdum))])
      disp(['size(y,1): ',num2str(ny)])
      error('NL2 constructor: size of h not consistent with nn')
   end
   try
      hdum=feval(h,[0 1],[x0 x0],zeros(ny,2),zeros(nu,2),th);
      if size(hdum,2)~=2
          error('NL2 constructor: size of h not consistent with nn')
      end
   catch
      disp('NL2 constructor warning: try to vectorize h for increased speed')
   end
end

if ~isnumericscalar(fs)
   error('NL2 constructor: fs must be a numeric scalar')
end
if fs<0 & ~isnan(fs)
   error('NL2 constructor: fs must be positive or NaN')
end
m.name=[];
m.desc=[];
ylabel=[]; ulabel=[]; xlabel=[]; thlabel=[];
for i=1:ny
   ylabel{i}=['y',num2str(i)];
end
for i=1:nx
   xlabel{i}=['x',num2str(i)];
end
for i=1:nu
   ulabel{i}=['u',num2str(i)];
end
for i=1:nth
   thlabel{i}=['th',num2str(i)];
end
tlabel='Time';
I=[];
P=[];
pe=[];
pv=[];
px0=[];
m.f=f; m.h=h; m.nn=nn; m.pv=pv; m.pe=pe;
m.px0=px0; m.x0=x0; m.fs=fs;
m.th=th; m.P=P; m.I=I;
m.xlabel=xlabel; m.thlabel=thlabel; m.ulabel=ulabel;
m.ylabel=ylabel; m.tlabel=tlabel;
end

function out=fieldread(X,arg)
%FIELDREAD allows reading the protected fields in the object
out=eval(arg);
end

function m=set.x0(m,x0);
   if ~isempty(x0) & length(x0)~=m.nn(1)
      error('NL2: cannot change dimension of x0')
   end
   m.x0=x0(:);
end

function m=set.th(m,th);
   if ~isempty(th) & length(th)~=m.nn(5)
      error('NL2: cannot change dimension of th')
   end
   m.th=th(:);
end

function m=set.P(m,P);
   if ~isempty(P) & length(P)~=m.nn(5)
      error('NL2: cannot change dimension of P')
   end
   if ~iscov(P)
      error('NL2: P is not a valid covariance')
   end
   m.P=P;
end

function m=set.I(m,I);
   if isempty(I)
      % nothing
   elseif  length(I)~=m.nn(1)+m.nn(5)
      error('NL2: cannot change dimension of I (square nx+nth matrix)')
   end
   if ~iscov(I)
      disp('NL2 warning: I is not a valid information matrix (symm pos semidefinit)')
   end
   m.I=I;
end

function m=set.pe(m,pe);
   if isempty(pe)
       m.pe=[];  % Remove the noise
   elseif isa(pe,'pdfclass')
       m.pe=pe;
       if length(pe)~=m.nn(4);  % length(e)=ny
          error('NL2 fieldwrite: length of symbolic pe must equal ny')
       end
   elseif iscov(pe)
       if size(pe,1)~=m.nn(4);  % length(e)=ny
          error('NL2 fieldwrite: size of covariance matrix pe must equal ny')
       end
       m.pe=ndist(zeros(m.nn(4),1),pe);
   else
      error('NL2 fieldwrite: pe must be either a symbolic distribution or a covariance matrix')
   end
end

function m=set.pv(m,pv);
   if isempty(pv)
       m.pv=[];  % Remove the noise
   elseif isa(pv,'pdfclass')
       m.pv=pv;
       if length(pv)~=m.nn(2);  % length(v)=nv
          error('NL2 fieldwrite: length of symbolic pv must equal nv')
       end
   elseif iscov(pv)
       if size(pv,1)~=m.nn(2);
          error('NL2 fieldwrite: size of covariance matrix pv must equal nv')
       end
       m.pv=ndist(zeros(m.nn(2),1),pv);
   else
      error('NL2 fieldwrite: pv must be either a symbolic distribution or a covariance matrix')
   end
end

function m=set.px0(m,px0);
   if isempty(px0)
       m.px0=[];  % Remove the noise
   elseif isa(px0,'pdfclass')
       m.px0=px0;
       if length(px0)~=m.nn(1);  % length(x0)=nx
          error('NL2 fieldwrite: length of symbolic px0 must equal nx')
       end
   elseif iscov(px0)
       if size(px0,1)~=m.nn(1);
          error('NL2 fieldwrite: size of covariance matrix px0 must equal nx')
       end
       m.px0=ndist(zeros(m.nn(1),1),px0);
   else
      error('NL2 fieldwrite: px0 must be either a symbolic distribution or a covariance matrix')
   end
end

function m=set.fs(m,fs);
   if isnumericscalar(fs) | fs<0
      m.fs=fs;
   else
      error('NL2 fieldwrite: fs must be a positive scalar or NaN')
   end
end

function m=set.xlabel(m,xlabel);
   if ~iscell(xlabel) & isstr(xlabel) & m.nn(1)==1
      m.xlabel={xlabel};   % Cellify
   elseif iscell(xlabel) & length(xlabel)==m.nn(1)
      m.xlabel=xlabel;
   elseif isempty(xlabel)
      m.xlabel=xlabel;
   else
      error(['NL2 fieldwrite: Field xlabel must be a cell of length nx, [], or string if nx=1'])
   end
end

function m=set.ylabel(m,ylabel);
   if ~iscell(ylabel) & isstr(ylabel) & m.nn(4)==1
      m.ylabel={ylabel};   % Cellify
   elseif iscell(ylabel) & length(ylabel)==m.nn(4)
      m.ylabel=ylabel;
   elseif isempty(ylabel)
      m.ylabel=ylabel;
   else
      error(['NL2 fieldwrite: Field ylabel must be a cell of length ny, [], or string if ny=1'])
   end
end

function m=set.ulabel(m,ulabel);
   if ~iscell(ulabel) & isstr(ulabel) & m.nn(3)==1
      m.ulabel={ulabel};   % Cellify
   elseif iscell(ulabel) & length(ulabel)==m.nn(3)
      m.ulabel=ulabel;
   elseif isempty(ulabel)
      m.ulabel=ulabel;
   else
      error(['NL2 fieldwrite: Field ulabel must be a cell of length nu, [], or string if nu=1'])
   end
end

function m=set.thlabel(m,thlabel);
   if ~iscell(thlabel) & isstr(thlabel) & m.nn(5)==1
      m.thlabel={thlabel};   % Cellify
   elseif iscell(thlabel) & length(thlabel)==m.nn(5)
      m.thlabel=thlabel;
   elseif isempty(thlabel)
      m.thlabel=thlabel;
   else
      error(['NL2 fieldwrite: Field thlabel must be a cell of length nth, [], or string if nth=1'])
   end
end

function m=set.tlabel(m,tlabel);
   if  isstr(tlabel)
      m.tlabel=tlabel;
   else
      error(['NL2 fieldwrite: Field tlabel must be a string'])
   end
end

function sys=arrayread(m,j,i)
%ARRAYREAD used to pick out sub-systems by indexing
%   mji=arrayread(m,j,i)
%   j is the row index/indices, corresponding to the outputs
%   i is the column index/indices, corresponding to the inputs
%   Example: s23=s(2,3) gives the SISO system from input 3 to output 2

if nargin<3
   error('NL2.ARRAYREAD: both input indeces in s(j,i) are required')
end
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
if isstr(i)
    if strcmp(i,':')
        i=1:nu;
    else
       error('NL2.ARRARYREAD: Inputs must be a vector or :')
    end
end
if isstr(j)
    if strcmp(j,':')
        j=1:ny;
    else
       error('NL2.ARRAYREAD: Inputs must be a vector or :')
    end
end
if any(i>nu)
     error('NL2.ARRAYREAD: Output (row) index larger than the number of outputs')
end
if any(j>ny)
     error('NL2.ARRAYREAD: Input (column) index larger than the number of inputs')
end
sysf=m.f;
sysh=m.h;
sysnn=m.nn;

if length(j)<ny % Delete outputs
   if ~isa(m.h,'inline')
      error('NL2.ARRAYREAD: Removing output is not possible for m-file definition of h')
   end
   I=eye(ny);
   sysh=inline([char(h),'*eye(j,:)']);
   sysnn(4)=length(j);
end
sys=nl2(sysf,sysh,sysnn);
sys.name=['sub-model of ',m.name];
if ~isempty(m.ylabel);
   sys.ylabel={m.ylabel{j}};
end
if ~isempty(m.ulabel);
   sys.ulabel={m.ulabel{i}};
end
if ~isempty(m.xlabel);
   sys.xlabel=m.xlabel;
end
end

function disp(m)
%DISP returns an ascii formatted version of the NL2 model
format='%11.2g';
if isa(m.f,'inline')
   fstr=formula(m.f);
else
   fstr=[m.f,'(t,x,u,th)'];
end
if isa(m.h,'inline')
   hstr=formula(m.h);
else
   hstr=[m.h,'(t,x,u,th)'];
end
if isa(m.pv,'pdfclass')
   pvstr=['  v = ',symbolic(m.pv)];
elseif isempty(m.pv)
   pvstr=[];
else
   pvstr=['   v = N(0,',num2str(m.pv,format),')'];
end
if length(pvstr)>20
   pvstr=[' + v'];
end
if isa(m.pe,'pdfclass')
   pestr=['  e = ',symbolic(m.pe)];
elseif isempty(m.pe)
   pestr=[];
else isnumeric(m.pe)
   pestr=['   e = N(0,',mat2strformat(m.pe,format),')'];
end
if length(pestr)>20
   pestr=[' + e'];
end
ystr='   y = ';
if isempty(m.px0) | all(all(cov(m.px0)==0))
   px0str=[];
else
   px0str=[' + N(0,',mat2strformat(cov(m.px0),format),')'];
end
if isnan(m.fs)
   xstr=' dx/dt = ';
else
   xstr='x[k+1] = ';
end

if isa(m,'sensormod')
   objname='SENSORMOD';
elseif isa(m,'sigmod')
   objname='SIGMOD';
else
   objname='NL2';
end
if isempty(m.name)
   disp([objname,' object'])
else
   disp([objname,' object: ',m.name])
end

if ~isempty(m.f) & ~strcmp(fstr,'x')
   ind=find(fstr==';');
   ind1=find(fstr=='[');
   ind2=find(fstr==']');
   n=length(ind);
   if n>0 & n==m.nn(1)-1 & ~strfind(fstr,'*x') % matrix with one row per measurement
      disp(['   x(1) = ',fstr(ind1(1)+1:ind(1)-1)])
      for i=2:n
         disp(['   x(',num2str(i),') = ',fstr(ind(i-1)+1:ind(i)-1)])
      end
      disp(['   x(',num2str(n+1),') = ',fstr(ind(n)+1:ind2(end)-1)])
   else
      fstr=[xstr,fstr,pvstr];
      try
         fstrx=inline2mat(m.f,xstr,pvstr);
         if size(fstrx,1)==m.nn(1)
            fstr=fstrx;
         end
      end
      [N,M]=size(fstr);
      k=0;
      nc=100;
      for k=1:floor(M/nc)
         disp(fstr(:,k*nc-nc+1:k*nc))
   %      disp(' ');
      end
      disp(fstr(:,floor(M/nc)*nc+1:end))
   end
end

if ~isempty(m.h)
   ind=find(hstr==';');
   ind1=find(hstr=='[');
   ind2=find(hstr==']');
   n=length(ind);
   if n>0 & n==m.nn(4)-1 & ~strfind(hstr,'*x') % matrix with one row per measurement
     disp(['   y(1) = ',hstr(ind1(1)+1:ind(1)-1)])
     for i=2:n
        disp(['   y(',num2str(i),') = ',hstr(ind(i-1)+1:ind(i)-1)])
     end
     disp(['   y(',num2str(n+1),') = ',hstr(ind(n)+1:ind2(end)-1)])
   else
      hstr=[ystr,hstr,pestr];
      try
         hstrx=inline2mat(m.h,ystr,pestr);
         if size(hstrx,1)==m.nn(4)
            hstr=hstrx;
         end
      end
      [N,M]=size(hstr);
      k=0;
      nc=100;
      for k=1:floor(M/nc)
         disp(hstr(:,k*nc-nc+1:k*nc))
      %   if m.nn(3)>1, disp(' '),end
      end
      disp(hstr(:,floor(M/nc)*nc+1:end))
   end
end
if ~isempty(m.f) & ~isa(m,'sigmod')
   disp(['   x0''  = ',mat2strformat(m.x0',format),'  ',px0str])
end
if ~isempty(m.P) & ~all(diag(m.P(m.nn(4)+1:end,m.nn(5)+1:end))==0)
    disp(['   std = [',num2str(sqrt(diag(m.P(m.nn(5)+1:end,m.nn(5)+1:end)))',format),']'])
end
if ~isempty(m.th)
    disp(['   th''  = ',mat2strformat(m.th',format)])
end
if ~isempty(m.P) & ~all(diag(m.P(1:m.nn(5),1:m.nn(5)))==0)
    disp(['   std = [',num2str(sqrt(diag(m.P(1:m.nn(5),1:m.nn(5))))',format),']'])
end
disp('')
if ~isempty(m.xlabel)
   str='  States:  ';
   for k=1:m.nn(1)
     str=[str, m.xlabel{k},'     '];
   end
   disp(str)
end
if ~isempty(m.ylabel)
   str='  Outputs: ';
   for k=1:m.nn(4)
     str=[str, m.ylabel{k},'     '];
   end
   disp(str)
end
if ~isempty(m.ulabel)
   str='  Inputs:  ';
   for k=1:m.nn(2)
     str=[str, m.ulabel{k},'     '];
   end
   disp(str)
end
if ~isempty(m.thlabel)
   str='  Param.: ';
   for k=1:m.nn(5)
     str=[str, m.thlabel{k},'     '];
   end
   disp(str)
end
end

function nn=size(m)
%SIZE returns the structure indeces nn=[nx,nu,nv,ny]
nn=m.nn;
end

function md=c2d(mc,fs)
%C2D discretizes a continuous time model to a discrete time model using Euler
%   md=c2d(mc,fs)
%   Euler sampling simply sets f_d(x)=x+1/fs*f_c(x) and pv_d=1/fs*pv_c

if isnan(mc.fs)
   fc=char(mc.f);
   nx=mc.nn(1);
   ind1=findstr(fc,'[');
   ind2=findstr(fc,';');
   ind=[ind1(1) ind2];
   if length(ind)~=nx, error('NL2.C2D string interpretation problems'), end
   ind=[ind length(fc)];
   fd=fc(1:ind(1));
   for k=1:nx
      fd=[fd,'x(',num2str(k),',:)+',num2str(1/fs),'*('];
      fd=[fd,fc(ind(k)+1:ind(k+1)-1),')',fc(ind(k+1))];
   end
   md=nl2(fd,mc.h,mc.nn);
   md.th=mc.th;
   md.pv=1/fs*mc.pv;   % Increase state noise
   md.pe=mc.pe;
   md.x0=mc.x0;
   md.px0=mc.px0;
   md.xlabel=mc.xlabel;
   md.ylabel=mc.ylabel;
   md.ulabel=mc.ulabel;
   md.name=mc.name;
   md.desc=mc.desc;
   md.fs=fs;
else
   error('NL2.C2D cannot change sampling interval of discrete time models')
end
end

function maug=augment(m,M)
%AUGMENT augments a motion model with M constant states
%   maug=augment(m,M)
%   Utility function for SLAM or recursive system identification problems
f=char(m.f);
f=['[',f];
for i=1:M;
   f=[f,'; x(',num2str(m.nn(1)+i),',:)'];
end
f=[f,']'];
maug=nl2(f,[],m.nn+[M 0 0 0]);
maug.fs=m.fs;
maug.x0=[m.x0;zeros(M,1)];
if ~isempty(m.px0)
   maug.px0=blkdiag(cov(m.px0),zeros(M));
end
if ~isempty(m.pv)
   maug.pv=blkdiag(cov(m.pv),zeros(M));
end
end

function ms=addsensor(m,s,varargin)
%ADDSENSOR adds (another) sensor to a NL2 model
%   ms=addsensor(m,s,Property1,Value1,...)
%
%   The nonl2inear model
%     x[k+1]=f(x[k])+v[k]
%       y[k]=h(x[k])+e[k]
%   where the sensor model might be empty,
%   is appended with another sensor model

%      y2[k]=h2(x[k])+e2[k]
%
%   Assumptions:
%      all parameters th are in s
%      the states in s.h are the first ones of the states in m.f
%      the sampling frequency in s is the same as in m
%      pe is Gaussian with covariance R
%
%   Property  Value  Description
%   ------------------------------------------
%
%   Example:
%      m=exmotion('ctcv2d');  % Motion model without sensor
%      s=exsensor('radar',1)  % Radar sensor
%      ms=addsensor(m,s);     % Motion model with one radar
%      mss=addsensor(ms,s);   % Motion model with two radars

nn=m.nn+[0 0 s.nn(3:4)];
if isa(m.h,'inline') & isa(s.h,'inline')
   h1=char(m.h);
   h2=char(s.h);
   ind1=find(h1==']');
   if isempty(ind1)
      h1=['[',h1];
   else
      h1=h1(1:ind1(end)-1);
   end
   ind2=find(h2=='[');
   if isempty(ind2)
      h2=[h2,']'];
   else
      h2=h2(ind2(1)+1:end);
   end
   h=[h1,';',h2];
   ms=nl2(m.f,h,nn);
   ms.pe=blkdiag(cov(m.pe),cov(s.pe));
   ms.ylabel={m.ylabel{:},s.ylabel{:}};
elseif isempty(m.h)
   ms=nl2(m.f,s.h,nn);
   if m.nn(5)>0 & s.nn(5)>0
     disp('Warning: both motion and sensor models contain parameters, check the result')
      ms.th={m.th{:},s.th{:}};
   elseif m.nn(5)>0 & s.nn(5)==0
      ms.th=m.th;
   elseif m.nn(5)==0 & s.nn(5)>0
      ms.th=s.th;
   end
   ms.pe=s.pe;
   ms.P=s.P;
   ms.ylabel=s.ylabel;
else
   error('nl2.addsensor: cannot mix mfiles and inline representations of h with this function')
end
ms.fs=m.fs;
ms.pv=m.pv;
ms.x0=m.x0;
ms.px0=m.px0;
ms.xlabel=m.xlabel;
ms.ulabel=m.ulabel;
ms.name=['Motion model: ',m.name,' Sensor model: ',s.name];
end

function ms=removesensor(m,ind)
%REMOVESENSOR removes the sensors ind in a NL2 model
%   ms=removesensor(m,ind)
%
%   In the nonl2inear model
%     x[k+1]=f(x[k])+v[k]
%       y[k]=h(x[k])+e[k]
%   the sensor on row ind in h and the corresponding measurements
%   y(ind) are removed
%
%   Assumptions:
%      pe is Gaussian with covariance R
%
%   Example:
%      m=exmotion('ctcv2d');  % Motion model without sensor
%      s=exsensor('radar',1)  % Radar sensor
%      ms1=addsensor(m,s);     % Motion model with one radar
%      ms2=removesensor(ms1,2); % Only range sensor kept

nn=m.nn-[0 0 length(ind) 0];
if nn(4)==0
   ms=nl2(m.f,[],nn);
   ms.fs=m.fs;
   ms.pv=m.pv;
   ms.x0=m.x0;
   ms.px0=m.px0;
   ms.xlabel=m.xlabel;
   ms.ulabel=m.ulabel;
elseif isa(m.h,'inline')
   h1=char(m.h);
   ind1=find(h1==';');
   ind2=find(h1=='[');
   if isempty(ind2)
      ind2=0;
   end
   ind3=find(h1==']')+1;
   if isempty(ind3)
      h1=[h1,';'];
      ind3=length(h1);
   else
      h1=[h1(1:end-1),';]'];  % add extra semi kolon at end
   end
   ind4=[ind2 ind1 ind3];
   for k=length(ind):-1:1;
       h1(ind4(ind(k))+1:ind4(ind(k)+1)-1)=[];
   end
   h1(end-1)=[];  % Remove last semi kolon
   ms=nl2(m.f,h1,nn);
   R=cov(m.pe);
   R(ind,:)=[];
   R(:,ind)=[];
   ms.pe=R;
   kind=1;
   for k=1:m.nn(4);
       if isempty(find(ind==k))
           ylab{kind}=m.ylabel{k};
           kind=kind+1;
       end
   end
   ms.ylabel=ylab;
   ms.fs=m.fs;
   ms.pv=m.pv;
   ms.x0=m.x0;
   ms.px0=m.px0;
   ms.xlabel=m.xlabel;
   ms.ulabel=m.ulabel;
   %ms.name=['Motion model: ',m.name,' Sensor model: ',s.name];
else
   error('nl2.removesensor: cannot mix mfiles and inline representations of h with this function')
end
end

function [x,V]=ekf(m,z,varargin)
%EKF implements the extended Kalman filter for state estimation
%   [x,V]=ekf(m,z,Property1,Value1,...)
%
%   m      NL2 object with model
%   z      SIG object with measurements
%   x      SIG object with state estimates
%          xhat=x.x and signal estimate yhat=x.y
%   V      Normalized sum of squared innovations, which should be a sequence
%          of chi2(nx) variables when the model is correct
%
%   User guidelines:
%   1. Increase state noise covariance Q to mitigate linearization errors in f
%   2. Increase noise covariance R to mitigate linearization errors in h
%
%   Property  Value  Description
%   ------------------------------------------
%   k          k>0 {0}   Prediction horizon:
%                        0 for filter (default)
%                        1 for one-step ahead predictor,
%   P0         {[]}      Initial covariance matrix
%                        Scalar value scales identity matrix
%                        Empty matrix gives a large identity matrix
%   x0         {[]}      Initial state matrix (overrides the value in m.x0)
%                        Empty matrix gives a zero vector
%   Q          {[]}      Process noise covariance (overrides the value in m.pv)
%                        Scalar value scales m.Q
%   R          {[]}      Measurement noise covariance (overrides the value in m)
%                        Scalar value scales m.pe
%
%   Example 1: non-linear model
%     m=exnl('ctpv2d'); % coordinated turn model
%     z=simulate(m,10); % ten seconds trajectory
%     zhat=ekf(m,z);    % EKF state estimation
%     xplot(z,zhat,'ind',[1 2],'conf',90);
%     figure, plot(z,zhat,'conf',90);
%
%   Example 2: linear model, comparing ekf on nl2 object with kalman on ss object
%     mss=exlti('cv2d');
%     mnl2=nl2(mss);
%     z=simulate(mss,10);
%     zhat1=kalman(mss,z);
%     zhat2=ekf(mnl2,z);
%     xplot2(z,zhat1,zhat2,'conf',90)


opt=struct('k',0,'P0',cov(m.px0),'x0',m.x0,'Q',cov(m.pv),'R',cov(m.pe));
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.EKF: only discrete time models')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
[N,nyz,nuz]=size(z);
if nyz~=ny
   error('NL2.EKF: Number of outputs ny must be the same in model and signal')
end
if nuz~=nu
   error('NL2.EKF: Number of inputs nu must be the same in model and signal')
end

y=z.y;
u=z.u;
t=z.t;

if isempty(m.pv)
   disp('NL2.EKF warning: m.pv is not defined, using v=0')
end
if isempty(m.pe)
   disp('NL2.EKF warning: m.pe is not defined, using e=0')
end
if isempty(opt.P0) | isnan(opt.P0)
   disp('NL2.EKF warning: px0 not defined, using a default value instead')
   P0=1e3*eye(nx);
elseif isnumericscalar(opt.P0)
    P0=opt.P0*eye(nx);
elseif iscov(opt.P0)
    P0=opt.P0;
elseif isa(opt.P0,'pdfclass')
    P0=cov(opt.P0);
else
   error('NL2.EKF: P0 must be empty, a scalar, a symmetric (nx,nx) matrix or a PDFCLASS object')
end

if isempty(opt.x0)  | isnan(opt.x0)
    % Use x0 from m.px0
    x0=m.x0;
elseif isnumeric(opt.x0) & size(opt.x0,1)==nx & size(opt.x0,2)==1
    x0=opt.x0;
else
   error('NL2.EKF: x0 must be empty or a (nx,1) vector')
end

if isempty(opt.Q)
    % Use Q from m.pv
    Q=cov(m.pv);
elseif isnumericscalar(opt.Q)
    Q=opt.Q*cov(m.pv);
elseif iscov(opt.Q)  & size(opt.Q,1)==nv & size(opt.Q,2)==nv
    Q=opt.Q;
else
   error('NL2.EKF: Q must be empty, a scalar, or a covariance (nx,nx) matrix')
end
if isempty(Q)
   Q=0;
end

if isempty(opt.R)
    % Use R from m
    R=cov(m.pe);
elseif isnumericscalar(opt.R)
    R=opt.R*cov(m.pe);
elseif iscov(opt.R) & size(opt.R,1)==ny & size(opt.R,2)==ny
    R=opt.R;
else
   error('NL2.EKF: Value of ''R'' must be empty, a scalar, or a covariance (ny,ny) matrix')
end
if isempty(R)
   error('NL2.EKF: Model must have a non-empty measurement noise covariance R')
end

if isnumericscalar(opt.k) & opt.k>=0
   %
else
   error('NL2.EKF: k must be a scalar k>=0')
end


xf=zeros(nx,N);
xp=zeros(nx,N+1);
Pf=zeros(N,nx,nx);
Pp=zeros(N+1,nx,nx);
yf=zeros(ny,N);
Pyf=zeros(N,ny,ny);
yp=zeros(ny,N);
Pyp=zeros(N,ny,ny);

xhat=x0;
P=P0;
if nu>0;
   u=z.u;
else
   u=zeros(N,0);
end
for k=1:N
   % Measurement update
   xlin=sig(zeros(1,ny),t(k),zeros(1,nu),xhat.');
   H=numgrad(m.h,2,t(k),xhat,zeros(ny,1),u(k,:),m.th);
   De=numgrad(m.h,3,t(k),xhat,zeros(ny,1),u(k,:),m.th);
   %H=nl2numgrad(m,xlin,'dhdx').'
   yphat=m.h(t(k),xhat,zeros(ny,1),u(k,:).',m.th);
   Rbar=De*R*De';
   Pyp(k,:,:)=H*P*H'+Rbar;
   Sinv=inv(H*P*H'+Rbar);
   K=P*H'*Sinv;
   P=P-K*H*P;
   P=0.5*(P+P');
   epsi=z.y(k,:).'-yphat;
   xhat=xhat+K*epsi;
   yfhat=m.h(t(k),xhat,zeros(ny,1),u(k,:),m.th);
   V(k)=epsi'*Sinv*epsi;

   xf(:,k)=xhat;
   Pf(k,:,:)=P;
   yf(:,k)=yfhat;
   yp(:,k)=yphat;
   Pyf(k,:,:)=H*P*H'+R;

   % Time update
   xlin=sig(zeros(1,ny),t(k),zeros(1,nu),xhat.');
   F=numgrad(m.f,2,t(k),xhat,zeros(nv,1),u(k,:).',m.th);
   G=numgrad(m.f,3,t(k),xhat,zeros(nv,1),u(k,:).',m.th);
   %F=nl2numgrad(m,xlin,'dfdx').';
   fxhat=feval(m.f,t(k),xhat,zeros(nv,1),u(k,:).',m.th);
   xhat=fxhat;
   P=F*P*F'+G*Q*G';
   xp(:,k+1)=xhat;
   Pp(k+1,:,:)=P;
end
xp(:,N+1)=[];
Pp(N+1,:,:)=[];
if opt.k==0
   x=sig(yf.',z.t,z.u,xf.',Pyf,Pf);
elseif opt.k==1
   x=sig(yp.',z.t,z.u,xp.',Pyp,Pp);
else
   error('NL2.EKF: k must be either 0 or 1')
end
x.fs=z.fs;
try
  x.xlabel=m.xlabel;
  x.ylabel=m.ylabel;
  x.ulabel=m.ulabel;
  x.name=m.name;
end
end

function I=fim(m,z,varargin)
%FIM computes the Fisher Information Matrix (FIM)
%   I=fim(m,z,Property1,Value1,...)
%
%   The FIM for y=h(x)+e is defined as
%       I=(dh/dx)^T*inv(R)*(dh/dx)
%   where the gradients are evaluated at the true state
%
%   m      NL2 object with model
%   z      State vector, or SIG object with true state in z.x
%   I      The (nx,nx) FIM  matrix
%
%   The syntax and options are identical to ekf.
%
%   Property  Value  Description
%   ------------------------------------------
%   R          {[]}      Measurement noise covariance (overrides the value in m)
%                        Scalar value scales m.R
%
%   Example:
%      m=nl2('x','[sqrt(x(1)^2+x(2)^2);atan2(x(2),x(1))]',[2 0 2 0],1);
%      m.pe=diag([0.1 0.01]);
%      for x1=1:3:10;
%        for x2=1:3:10;
%           I=fim(m,[x1;x2]);
%           mm=ndist([x1;x2],inv(I));
%           plot2(mm,'legend','','levels',1), hold on
%        end
%      end

opt=struct('R',cov(m.pe));
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.FIM: only discrete time models')
end
nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
if isa(z,'sig')
    if isempty(z.x)
        error('NL2.FIM: the input SIG object z must have a state x')
    end
    x=z.x;
    t=z.t;
    [N,nyz,nuz]=size(z);
else
    x=z;
    t=1;
    [N,nxdata]=size(x);
    if nx~=nxdata
       if N==nx
          x=x';
          N=1;
       else
           error(['NL2.FIM: the input state vector z must must have ',num2str(nx),' elements'])
       end
    end
end

if N>1
    error('NL2.FIM: FIM is only evaluated for one state at the time')
end

if isempty(opt.R)
    % Use R from m
    R=cov(m.pe);
elseif isnumericscalar(opt.R)
    R=opt.R*cov(m.pe);
elseif iscov(opt.R) & size(opt.R,1)==ny & size(opt.R,2)==ny
    R=opt.R;
else
   error('NL2.FIM: Value of ''R'' must be empty, a scalar, or a covariance (ny,ny) matrix')
end
if isempty(R)
   error('NL2.FIM: Model must have a non-empty measurement noise covariance R')
end

xsig=sig(zeros(1,ny),t,zeros(1,nu),x);
H=nl2numgrad(m,xsig,'dhdx').';
I=H'*inv(R)*H;
end


function [p,t,lambda]=pd(m,z,pfa,varargin)
%PD computes the probability of detection (PD) using GLRT
%   [p,t,lambda]=pd(m,z,pfa,Property1,Value1,...)
%
%   Consider the hypothesis test
%      H0: y = e
%      H0: y = h(x)+e
%   The generalized likelihood ratio test (GLRT) gives the test statistic
%       T(y) = chi2(nx) under H0
%       T(y) = ncchi2(nx,lambda) under H1
%   where
%       lambda = x'*J(x)*x
%   and J(x) is the Fisher information matrix (FIM).
%   The false alarm probabilty (PFA) is given by PFA=Prob(T(y)>t|H0),
%   which implicitely defines the test threshold t.
%   The probabilty of detection (PD) is given PD=Prob(T(y)>t|H1),
%
%   m      NL2 object with model
%   z      State vector, or SIG object with true state in z.x
%   I      The (nx,nx) PD  matrix
%
%   The syntax and options are identical to fim.
%
%   Property  Value  Description
%   ------------------------------------------
%   R          {[]}      Measurement noise covariance (overrides the value in m)
%                        Scalar value scales m.R
%
%   Example:
%   % PD, t, lambda for the test H0: y=e, H1: y=x+e, with PFA=0.01
%     nx=2;
%     m1=nl2('x','x',[nx 0 nx 0],1);
%     m1.pe=eye(nx);
%     [p(nx),t(nx),lambda(nx)]=pd(m1,1*ones(nx,1),0.01);
%   % PD(X,Y) for radar observations
%      m=nl2('x','[sqrt(x(1)^2+x(2)^2);atan2(x(2),x(1))]',[2 0 2 0],1);
%      m.pe=diag([0.1 0.01]);
%      x1=0.5:0.5:3;
%      x2=0.5:0.5:3;
%      for i=1:length(x1)
%         for j=1:length(x2)
%            p(i,j)=pd(m,[x1(i);x2(j)]);
%        end
%      end
%      p
%      contour(x1,x2,p)

if nargin<2
    error('NL2.PD: at least two input arguments required')
end
if nargin<3
   pfa=0.01;
end
if ~isnumericscalar(pfa) | pfa<=0 | pfa>=1
    error('NL2.PD: 0<pfa<1')
end

if isnan(m.fs)
    error('NL2.PD: only discrete time models')
end
nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
if isa(z,'sig')
    if isempty(z.x)
        error('NL2.PD: the input SIG object z must have a state x')
    end
    x=z.x;
    t=z.t;
    [N,nyz,nuz]=size(z);
else
    x=z;
    t=1;
    [N,nxdata]=size(x);
    if nx~=nxdata
       if N==nx
          x=x';
          N=1;
       else
           error(['NL2.PD: the input state vector z must must have ',num2str(nx),' elements'])
       end
    end
end

if N>1
    error('NL2.PD: PD is only evaluated for one state at the time')
end

nx=m.nn(1);
t=erfinv(chi2dist(nx),1-pfa);
I=fim(m,z,varargin{:});
lambda=x*I*x';
p=1-erf(ncchi2dist(nx,lambda),t);
end



function x=crlb(m,z,varargin)
%CRLB computes the parametric Cramer Rao Lower Bound for state estimation
%   x=crlb(m,z,Property1,Value1,...)
%
%   m      NL2 object with model
%   z      SIG object with true states in z.x
%   x      SIG object with x.x=z.x and x.y=h(z.x), where the covariances are
%          replaced by the corresponding Cramer-Rao lower bound on x and y
%
%   The syntax and options are identical to ekf.
%
%   Property  Value  Description
%   ------------------------------------------
%   k          k>0 {0}   Prediction horizon:
%                        0 for filter (default)
%                        1 for one-step ahead predictor,
%   P0         {[]}      Initial covariance matrix
%                        Scalar value scales identity matrix
%                        Empty matrix gives a large identity matrix
%   x0         {[]}      Initial state matrix (overrides the value in m.x0)
%                        Empty matrix gives a zero vector
%   Q          {[]}      Process noise covariance (overrides the value in m.Q)
%                        Scalar value scales m.Q
%   R          {[]}      Measurement noise covariance (overrides the value in m)
%                        Scalar value scales m.R
%
%   Example 1: non-linear model
%     m=exnl('ctpv2d'); % coordinated turn model
%     z=simulate(m,10); % ten seconds trajectory
%     zekf=ekf(m,z);
%     zcrlb=crlb(m,z);
%     xplot2(z,zekf,zcrlb,'conf',90)
%     figure, plot(z,zekf,zcrlb,'conf',90);
%
%   Example 2: linear model, comparing ekf on nl2 object with kalman on ss object
%     mss=exlti('cv2d');
%     mnl2=nl2(mss);
%     z=simulate(mss,10);
%     zhat1=kalman(mss,z);
%     zhat2=ekf(mnl2,z);
%     zcrlb=ekf(mnl2,z);
%     xplot2(z,zhat1,zhat2,zcrlb,'conf',90)

opt=struct('k',0,'P0',cov(m.px0),'x0',mean(m.px0),'Q',cov(m.pv),'R',cov(m.pe));
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.CRLB: only discrete time models')
end
if isempty(z.x)
    error('NL2.CRLB: the input SIG object z must have a state x')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
[N,nyz,nuz]=size(z);
%if nyz~=ny
%   error('NL2.CRLB: Number of outputs ny must be the same in model and signal')
%end
%if nuz~=nu
%   error('NL2.CRLB: Number of inputs nu must be the same in model and signal')
%end

y=z.y;
u=z.u;
t=z.t;

if isempty(opt.P0) | isnan(opt.P0)
    P0=1e3*eye(nx);
elseif isnumericscalar(opt.P0)
    P0=opt.P0*eye(nx);
elseif iscov(opt.P0)
    P0=opt.P0;
elseif isa(opt.P0,'pdfclass')
    P0=cov(opt.P0);
else
   error('NL2.CRLB: P0 must be empty, a scalar, a symmetric (nx,nx) matrix or a PDFCLASS object')
end

if isempty(opt.x0) | isnan(opt.x0)
    % Use x0 from m.px0
    x0=m.x0;
elseif isnumeric(opt.x0) & size(opt.x0,1)==nx & size(opt.x0,2)==1
    x0=opt.x0;
else
   error('NL2.CRLB: x0 must be empty or a (nx,1) vector')
end

if isempty(opt.Q)  | isnan(opt.Q)
    % Use Q from m.pv
    Q=cov(m.pv);
elseif isnumericscalar(opt.Q)
    Q=opt.Q*cov(m.pv);
elseif iscov(opt.Q)  & size(opt.Q,1)==nv & size(opt.Q,2)==nv %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q=opt.Q;
else
   error('NL2.CRLB: Q must be empty, a scalar, or a covariance (nx,nx) matrix')
end
if isempty(Q) | isnan(Q)
   Q=0;
end

if isempty(opt.R)
    % Use R from m
    R=cov(m.pe);
elseif isnumericscalar(opt.R)
    R=opt.R*cov(m.pe);
elseif iscov(opt.R) & size(opt.R,1)==ny & size(opt.R,2)==ny
    R=opt.R;
else
   error('NL2.CRLB: Value of ''R'' must be empty, a scalar, or a covariance (ny,ny) matrix')
end
if isempty(R)
   error('NL2.CRLB: Model must have a non-empty measurement noise covariance R')
end
if all(diag(R)==0)
   error('NL2.CRLB: Model must have a positive definite measurement noise covariance R')
end

if isnumericscalar(opt.k) & opt.k>=0
   %
else
   error('NL2.CRLB: k must be a scalar k>=0')
end

xf=zeros(nx,N);
xp=zeros(nx,N);
Pf=zeros(N,nx,nx);
Pp=zeros(N,nx,nx);
yf=zeros(ny,N);
Pyf=zeros(N,ny,ny);
yp=zeros(ny,N);
Pyp=zeros(N,ny,ny);

xhat=m.x0;
P=P0;
for k=1:N
   if nu>0
      uk=u(k,:).';
   else
      uk=[];
   end
   % Measurement update
   xhat=z.x(k,:).';  % Replace xhat with true state
   yfhat=feval(m.h,t(k),xhat,uk,m.th);
   xlin=sig(zeros(1,ny),t(k),zeros(1,nu),xhat.');
   C=nl2numgrad(m,xlin,'dhdx').';
   K=P*C'*inv(C*P*C'+R);
   P=P-K*C*P;
   P=0.5*(P+P');

   xf(:,k)=xhat;
   Pf(k,:,:)=P;
   yf(:,k)=yfhat;
   Pyf(k,:,:)=C*P*C'+R;

   % Time update
   xlin=sig(zeros(1,ny),t(k),zeros(1,nu),xhat.');
   A=nl2numgrad(m,xlin,'dfdx').';
   P=A*P*A'+Q;
   xp(:,k)=xhat;
   yp(:,k)=yfhat;
   Pp(k,:,:)=P;
   Pyp(k,:,:)=C*P*C'+R;
end
if opt.k==0
   x=sig(yf.',z.t,z.u,xf.',Pyf,Pf);
elseif opt.k==1
   x=sig(yp.',z.t,z.u,xp.',Pyp,Pp);
else
   error('NL2.CRLB: k must be either 0 or 1')
end
x.fs=z.fs;
try
  x.xlabel=m.xlabel;
  x.ylabel=m.ylabel;
  x.ulabel=m.ulabel;
  x.name=m.name;
end
end


function [x,V]=ukf(varargin)
%UKF implements the unscented Kalman filter
%   [x,V]=ukf(m,z,Property1,Value1,...)
%   The function calls nl2tf, so consult its help text
[x,V]=nl2tf(varargin{:},'tup','uteval','mup','uteval');
end

function [x,V]=ekf2(varargin)
%EKF2 implements the second order extended Kalman filter
%   [x,V]=ekf2(m,z,Property1,Value1,...)
%   The function calls nl2tf, so consult its help text
[x,V]=nl2tf(varargin{:},'tup','tt2eval','mup','tt2eval');
end

function [x,V]=ekf1(varargin)
%EKF1 implements the second order extended Kalman filter
%   [x,V]=ekf1(m,z,Property1,Value1,...)
%   ekf1 is another implementation than ekf, but should give the same result.
%   The function calls nl2tf, so consult its help text
[x,V]=nl2tf(varargin{:},'tup','tt1eval','mup','tt1eval');
end

function [x,V]=nl2tf(m,z,varargin)
%NL2TF implements the UKF and EKF without Ricatti equations
%   [x,V]=nl2tf(m,z,Property1,Value1,...)
%
%   m      NL2 object with model
%   z      SIG object with measurements
%   x      SIG object with state estimates
%          xhat=x.x and signal estimate yhat=x.y
%   V      Normalized squared innovation ideally chi2(nx)
%
%   The algorithm works as follows:
%   1. Time update:
%      a. Let xbar = [x;v] = N([xhat;0];[P,0;0,Q])
%      b. Transform approximation of x(k+1) = f(x,u,v)
%         gives xhat,P
%   2. Measurement update:
%      a. Let xbar = [x;e] = N([xhat;0];[P,0;0,R])
%      b. Transform approximation of z(k) = [x;y] = [x;h(x,u,v)]
%         provides zhat=[xhat;yhat] and Pz=[Pxx Pxy;Pyx Pyy]
%      c. The Kalman gain is K=Pxy*inv(Pyy)
%      d. xhat = xhat+K*(y-yhat)
%         P = P-K*Pyy*K'
%
%   The transform in 1b and 2b can be chosen arbitrarily as
%   uteval, tt1eval, tt2eval and mceval in the ndist object
%
%   Note: the NL2 object must be a function of indexed states, so always
%   write for instance x(1,:) or x(1:end,:), even for scalar systems.
%   The reason is that the state vector is augmented, so any unindexed x
%   will cause errors.
%
%   User guidelines:
%   1. Increase state noise covariance Q to mitigate linearization errors in f
%   2. Increase noise covariance R to mitigate linearization errors in h
%   3. Avoid very large values of P0 and Q (which can be used for KF and EKF)
%
%   Property   Value      Description
%   ------------------------------------------
%   k          k>0 {0}    Prediction horizon:
%                         0 for filter (default)
%                         1 for one-step ahead predictor,
%   P0         {[]}       Initial covariance matrix
%                         Scalar value scales identity matrix
%                         Empty matrix gives a large identity matrix
%   x0         {[]}       Initial state matrix (overrides the value in m.x0)
%                         Empty matrix gives a zero vector
%   Q          {[]}       Process noise covariance (overrides m.Q)
%                         Scalar value scales m.Q
%   R          {[]}       Measurement noise covariance (overrides m.R)
%                         Scalar value scales m.R
%   tup        {'uteval'} The unscented Kalman filter (UKF)
%              'tt1eval'  The extended Kalman filter (EKF)
%              'tt2eval'  The second order extended Kalman filter
%              'mceval'   The Monte Carlo KF
%   mup        {'uteval'} The unscented Kalman filter (UKF)
%              'tt1eval'  The extended Kalman filter (EKF)
%              'tt2eval'  The second order extended Kalman filter
%              'mceval'   The Monte Carlo KF
%  ukftype     'ut1',{'ut2'}|'ct'  Standard, modified UT, or cubature transform
%  ukfpar      {[]}       Parameters in UKF
%                         For ut1, par=w0 {w0=1-n/3}
%                         For ut2, par=[beta,alpha,kappa] {[2 1e-3 0]}
%                         For ct,  par=[a] with default [1]
%  NMC          {100}     Number of Monte Carlo samples for mceval
%
%   Example 1: non-linear model
%     m=exnl('ctpv2d');  % coordinated turn model
%     z=simulate(m,10);  % ten seconds trajectory
%     zukf=nl2tf(m,z); % UKF state estimation
%     zekf=nl2tf(m,z,'tup','tt1','mup','tt1'); % EKF variant
%     xplot(z,zukf,zekf,'ind',[1 2],'conf',90);


opt=struct('k',0,'P0',cov(m.px0),'x0',m.x0,'Q',cov(m.pv),'R',cov(m.pe),'tup','uteval','mup','uteval','ukftype','ut2','ukfpar',[],'NMC',100);
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.NL2TF: only discrete time models')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
[N,nyz,nuz]=size(z);
if nyz~=ny
   error('NL2.NL2TF: Number of outputs ny must be the same in model and signal')
end
if nuz~=nu
   error('NL2.NL2TF: Number of inputs nu must be the same in model and signal')
end

y=z.y;
u=z.u;
t=z.t;


if isempty(m.pv)
   disp('NL2.NL2TF warning: m.pv is not defined, using v=0')
end
if isempty(m.pe)
   disp('NL2.NL2TF warning: m.pe is not defined, using e=0')
end
if isempty(opt.P0)  | isnan(opt.P0)
    % Use P0 from m.px0
    disp('NL2.NL2TF warning: px0 not defined, using a default value instead')
    P0=1e3*eye(nx);
elseif isnumericscalar(opt.P0)
    P0=opt.P0*eye(nx);
elseif iscov(opt.P0)
    P0=opt.P0;
elseif isa(opt.P0,'pdfclass')
    P0=cov(opt.P0);
else
   error('NL2.NL2TF: P0 must be empty, a scalar, a symmetric (nx,nx) matrix or a PDFCLASS object')
end

if isempty(opt.x0) | isnan(opt.x0)
    % Use x0 from m.px0
    x0=m.x0;
elseif isnumeric(opt.x0) & size(opt.x0,1)==nx & size(opt.x0,2)==1
    x0=opt.x0;
else
   error('NL2.NL2TF: x0 must be empty or a (nx,1) vector')
end

if isempty(opt.Q)
    % Use Q from m.pv
    Q=cov(m.pv);
elseif isnumericscalar(opt.Q)
    Q=opt.Q*cov(m.pv);
elseif iscov(opt.Q)  & size(opt.Q,1)==nx & size(opt.Q,2)==nx
    Q=opt.Q;
else
   error('NL2.NL2TF: Q must be empty, a scalar, or a covariance (nx,nx) matrix')
end
if isempty(Q)
   Q=0;
end

if isempty(opt.R)
    % Use R from m
    R=cov(m.pe);
elseif isnumericscalar(opt.R)
    R=opt.R*cov(m.pe);
elseif iscov(opt.R) & size(opt.R,1)==ny & size(opt.R,2)==ny
    R=opt.R;
else
   error('NL2.NL2TF: Value of ''R'' must be empty, a scalar, or a covariance (ny,ny) matrix')
end
if isempty(R)
   error('NL2.NL2TF: Model must have a non-empty measurement noise covariance R')
end

if isnumericscalar(opt.k) & opt.k>=0
   %
else
   error('NL2.NL2TF: k must be a scalar k>=0')
end

if strncmpi(opt.mup,'tt1',3) | strcmpi(opt.mup,'taylor1')
   optpar=cell(0);
   mup='tt1eval';
elseif strncmpi(opt.mup,'tt2',3) | strcmpi(opt.mup,'taylor2')
   optpar=cell(0);
   mup='tt2eval';
elseif strncmpi(opt.mup,'ut',2)
   optpar=cell(1,2);
   optpar{1}=opt.ukftype;
   optpar{2}=opt.ukfpar;
   mup='uteval';
elseif strncmpi(opt.mup,'mc',2)
   optpar=cell(1);
   optpar{1}=opt.NMC;
   mup='mceval';
else
   error(['Unknown option for mup: ',opt.mup])
end

if strncmpi(opt.tup,'tt1',3) | strcmpi(opt.tup,'taylor1')
   optpar=cell(0);
   tup='tt1eval';
elseif strncmpi(opt.tup,'tt2',3) | strcmpi(opt.tup,'taylor2')
   optpar=cell(0);
   tup='tt2eval';
elseif strncmpi(opt.tup,'ut',2)
   optpar=cell(1,2);
   optpar{1}=opt.ukftype;
   optpar{2}=opt.ukfpar;
   tup='uteval';
elseif strncmpi(opt.tup,'mc',2)
   optpar=cell(1);
   optpar{1}=opt.NMC;
   tup='mceval';
else
   error(['Unknown option for tup: ',opt.tup])
end


xf=zeros(nx,N);
xp=zeros(nx,N);
Pf=zeros(N,nx,nx);
Pp=zeros(N,nx,nx);
yf=zeros(ny,N);
Pyf=zeros(N,ny,ny);
yp=zeros(ny,N);
Pyp=zeros(N,ny,ny);
if nu>0;
   u=z.u;
else
   u=zeros(N,0);
end
y=z.y.';
t=z.t.';
xhat=x0;
P=P0;
for k=1:N
   % Measurement update
   Xbar=ndist([xhat;zeros(ny,1)],[P zeros(nx,ny);zeros(ny,nx) R]);
   h=char(m.h);
   fbar=['[x(1:',num2str(nx),',:);',h,'+x(',num2str(nx),'+1:end,:)]'];
   fbar=inline(fbar,'x','t','u','th');
   Z=feval(mup,Xbar,fbar,optpar{:},t(k),u(k,:).',m.th);
   Ez=E(Z);
   Pz=cov(Z);
   yfhat=y(:,k);
   yphat=Ez(nx+1:end);
   P=Pz(1:nx,1:nx);
   Pyy=Pz(nx+1:end,nx+1:end);
   Pyyinv=inv(Pyy);
   K=Pz(1:nx,nx+1:end)*Pyyinv;
   epsi=y(:,k)-yphat;
   xhat=xhat+K*epsi;
   P=P-K*Pyy*K';
   P=0.5*(P+P');
   V(k)=epsi'*Pyyinv*epsi;

   xf(:,k)=xhat;
   Pf(k,:,:)=P;
   yf(:,k)=yfhat;
   yp(:,k)=yphat;
   Pyp(k,:,:)=Pz(nx+1:end,nx+1:end);
   Pyf(k,:,:)=Pz(nx+1:end,nx+1:end);  % XXX incorrect

   % Time update
   [Uq,Dq]=svd(Q);
   nv=rank(Q);
   Qbar=Dq(1:nv,1:nv);
   Bv=Uq(:,1:nv);
   Xbar=ndist([xhat;zeros(nv,1)],[P zeros(nx,nv);zeros(nv,nx) Qbar]);
   h=char(m.f);
   fbar=[h,'+',mat2strformat(Bv),'*x(',num2str(nx),'+1:end,:)'];
   fbar=inline(fbar,'x','t','u','th');
   Z=feval(tup,Xbar,fbar,optpar{:},t(k),u(k,:).',m.th);
   xhat=E(Z);
   P=cov(Z);
   xp(:,k)=xhat;
   Pp(k,:,:)=P;
end
if opt.k==0
   x=sig(yf.',t,u,xf.',Pyf,Pf);
elseif opt.k==1
   x=sig(yp.',t,u,xp.',Pyp,Pp);
else
   error('NL2.NL2TF: k must be either 0 or 1')
end
x.fs=z.fs;
try
  x.xlabel=m.xlabel;
  x.ylabel=m.ylabel;
  x.ulabel=m.ulabel;
  x.name=m.name;
end
end


function zhat=pf(m,z,varargin)
%PF implements the particle filter for state estimation
%   zhat=pf(m,z,Property1,Value1,...)
%
%   m      NL2 object with model
%   z      SIG object with measurements
%   x      SIG object with state estimates
%          xhat=x.x and signal estimate yhat=x.y
%
%   Property   Value      Description
%   ------------------------------------------
%   Np         Np>0 {100} Number of particles
%   k          k=0,1      Prediction horizon:
%                         0 for filter (default)
%                         1 for one-step ahead predictor,
%   proposal  {'sir'}     SIR, standard bootstrap PF (for low SNR)
%             'opt'       Approximation of optimal proposal (for high SNR)
%   sampling  {'simple'}  Standard multinomial sampling
%            'systematic'
%            'residual'
%            'stratified'
%   animate  {'off'},'on' Animation of particles on or off
%   ind      {[]},ind     Animate states x(ind), [] means all states
%   delay    T>=0         pause(T) after each plot
%   type     {'xplot'}    Plot the states x(ind) as subplots of stem
%            'plot'       Plot the outputs y(ind) as subplots of stem
%            'xplot2'     Plot of states x(ind(1)) to x(ind(2))
%
%   Example 1: non-linear model
%     m=exnl('ctpv2d'); % Coordinated turn model
%     z=simulate(m,10); % Ten seconds trajectory
%     zhat=pf(m,z);     % PF state estimation
%     xplot(z,zhat,'ind',[1 2],'conf',90);
%     figure, plot(z,zhat,'conf',90);
%
%   Example 2: linear model, comparing ekf on nl2 object with kalman on ss object
%     mss=exlti('cv2d');
%     mnl2=nl2(mss);
%     z=simulate(mss,10);
%     zhat1=kalman(mss,z);
%     zhat2=pf(mnl2,z);
%     xplot2(z,zhat1,zhat2,'conf',90)

opt=struct('k',0,'Np',100,'sampling','simple','proposal','sir','animate','off','ind',[],'delay',0,'type','xplot','Xlim',[],'Ylim',[]);
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.PF: only discrete time models')
end

if ~isnumericscalar(opt.k) | opt.k<0 | opt.k>1
   error('NL2.PF: k must be 0 or 1')
end
if ~isnumericscalar(opt.Np) | opt.Np<1
   error('NL2.PF: Np must a positive scalar')
end

if isempty(opt.ind)
   ind=':';
elseif ~isnumeric(opt.ind)
   error('NL2.PF: value for argument animate must be numeric')
else
   ind=opt.ind;
end
if ~isnumeric(opt.delay) | opt.delay<0
   error('NL2.PF: delay must be numeric and positive')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
[N,nyz,nuz]=size(z);
Np=opt.Np;
if nyz~=ny
   error('NL2.PF: Number of outputs ny must be the same in model and signal')
end
if nuz~=nu
   error('NL2.PF: Number of inputs nu must be the same in model and signal')
end

y=z.y.';
if nu>0;
   u=z.u.';
else
   u=zeros(0,N);
end
t=z.t;
xhat=zeros(N,nx);
Px=zeros(Np,N,nx);
yhat=zeros(N,ny);
Py=zeros(Np,N,ny);
w=1/Np*ones(Np,1);
xp=ones(Np,1)*m.x0.';
if isa(m.px0,'pdfclass')
   xp=xp + rand(m.px0,Np);
end

for k=1:N;
   % Prediction y(k|k-1)
      try  % vectorized call
         yp=feval(m.h,t(k),xp.',u(:,k).',m.th).';
      catch
         for i=1:Np
            yp(i,:)=feval(m.h,t(k),xp(i,:).',u(:,k),m.th).';
         end
      end
   if opt.k==1
      xhat(k,:)=mean(xp);
      Px(:,k,:)=xp;
      yhat(k,:)=mean(yp);
      Py(:,k,:)=yp;
   end
   % Measurement update of x(k|k)
   if strncmpi(opt.proposal,'sir',3)
      % SIR (Bootstrap) PF
      w=w.*pdf(m.pe,repmat(y(:,k).',Np,1)-yp);     % Likelihood
   elseif strncmpi(opt.proposal,'opt',3)
      % Approximation of the optimal proposal
      Q=cov(m.pv);
      R=cov(m.pe);
      xk=sig(xhat(k,:),t(k),u(:,k).',xhat(k,:));
      H=nl2numgrad(m,xk,'dhdx').';
      epsi=repmat(y(:,k).',Np,1)-yp;     % prediction error
      w=w.*pdf(ndist(zeros(ny,1),R+H*Q*H'),epsi);     % Likelihood approximation
   end
   % Animate
   if strcmp(opt.animate,'on')
      if strcmp(opt.type,'xplot')
         if strcmp(ind,':'), ind=1:nx; end
         for i=1:length(ind)
            subplot(length(ind),1,i)
            %hist(xp(:,ind(i)))
            stem(xp(:,ind(i)),w(:))
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            if ~isempty(z.x) % True value available
               hold on
               ax=axis;
               plot(z.x(k,ind(i))*[1 1],[ax(3);ax(4)],'r')
               hold off
               title([z.xlabel{ind(i)},' at time ',num2str(t(k))])
            end
         end
         drawnow
         pause(opt.delay)
      elseif strcmp(opt.type,'xplot2')
         if strcmp(ind,':'), ind=[1 2]; end
         if length(ind)~=2
            error('NL2.PF: ind must have length 2 in xplot2')
         end
            hp=plot(xp(:,ind(1)),xp(:,ind(2)),'b.');
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            if ~isempty(opt.Ylim), set(gca,'Ylim',opt.Ylim), end
            if ~isempty(z.x) % True value available
               hold on
               plot(z.x(k,ind(1)),z.x(k,ind(2)),'*r')
               hold off
               title(['Time ',num2str(t(k))])
               xlabel(z.xlabel{ind(1)})
               ylabel(z.xlabel{ind(2)})
            end
         drawnow
         pause(opt.delay)
         if k<N
            set(hp,'Xdata',zeros(1,Np),'Ydata',zeros(1,Np))
         end
      elseif strcmp(opt.type,'plot')
         if strcmp(ind,':'), ind=1:ny; end
         for i=1:length(ind)
            subplot(length(ind),1,i)
            %hist(xp(:,ind(i)))
            stem(yp(:,ind(i)),w(:))
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            hold on
            ax=axis;
            plot(z.y(k,ind(i))*[1 1],[ax(3);ax(4)],'r')
            hold off
            title([z.ylabel{ind(i)},' at time ',num2str(t(k))])
         end
         drawnow
         pause(opt.delay)
      end
   end
   [xp,w]=resampling(xp,w,opt.sampling);  % Resampling
   if opt.k==0
      % Compute y(k|k)
        try  % vectorized call
           yp=feval(m.h,t(k),xp.',u(:,k).',m.th).';
        catch
           for i=1:Np
              yp(i,:)=feval(m.h,t(k),xp(i,:).',u(:,k),m.th).';
           end
        end
      xhat(k,:)=mean(xp);
      Px(:,k,:)=xp;
      yhat(k,:)=mean(yp);
      Py(:,k,:)=yp;
   end
   % Time update of x(k+1|k)
   try  % vectorized call
      xp=feval(m.f,t(k),xp.',u(:,k),m.th).';
   catch
      for i=1:Np
         xp(i,:)=feval(m.f,t(k),xp(i,:).',u(:,k),m.th).';
      end
   end
   if strncmpi(opt.proposal,'sir',3)
      % SIR (Bootstrap) PF
      v=rand(m.pv,Np);
   elseif strncmpi(opt.proposal,'opt',3)
      % Approximation of the optimal proposal
      if k<N;
         Q=cov(m.pv);
         R=cov(m.pe);
         xk=sig(xhat(k,:),t(k+1),u(:,k+1).',xhat(k,:));
         H=nl2numgrad(m,xk,'dhdx').';
         Qbar=pinv(H'*pinv(R)*H+pinv(Q));
         try  % vectorized call
            yp=feval(m.h,t(k+1),xp.',u(:,k+1).',m.th).';
         catch
            for i=1:Np
               yp(i,:)=feval(m.h,t(k),xp(i,:).',u(:,k),m.th).';
            end
         end
         epsi=repmat(y(:,k).',Np,1)-yp;     % prediction error
         K=Q*H'*inv(H*Q*H'+R);
         Xprop=ndist(zeros(nx,1),Qbar);
         v2=rand(Xprop,Np);  % Same distribution for each particle
         v1=(K*epsi.').';    % Different mean for each particle
         v=v1+v2;
         w=w.*pdf(m.pv,v)./pdf(Xprop,v);
         %         epsi,xp,pause
      else
         v=zeros(size(xp)); % Do nothing last time update
      end
   end
   xp=xp+v;
end
%plot(xp(1,:),xp(2,:),'r.')
zhat=sig(yhat,z.t,u.',xhat,Py,Px);
zhat.fs=z.fs;
try
  zhat.xlabel=m.xlabel;
  zhat.ylabel=m.ylabel;
  zhat.ulabel=m.ulabel;
  zhat.name=m.name;
end
end


function zhat=pmf(m,z,xrange,varargin)
%PMF implements the point mass filter for state estimation
%   zhat=pmf(m,z,xrange,Property1,Value1,...)
%
%   m      NL2 object with model
%   z      SIG object with measurements
%   x      SIG object with state estimates
%          xhat=x.x and signal estimate yhat=x.y
%   xrange (2,nx) matrix with min and max values for each x dimension
%          The grid consists of Np uniformly spaced grid points
%
%   Property   Value      Description
%   ------------------------------------------
%   Np         Np>0 {100} Number of particles
%   k          k=0,1      Prediction horizon:
%                         0 for filter (default)
%                         1 for one-step ahead predictor,
%   animate  {'off'},'on' Animation of particles on or off
%   ind      {[]},ind     Animate states x(ind), [] means all states
%   delay    T>=0         pause(T) after each plot
%   type     {'xplot'}    Plot the states x(ind) as subplots of stem
%            'plot'       Plot the outputs y(ind) as subplots of stem
%            'xplot2'     Plot of states x(ind(1)) to x(ind(2))
%
%   Example: random walk model, radar
%     m=exnl('cp2dradar');
%     xs=[0:9;0:2:18];  ys=m.h(0,xs,[],[]);
%     y=sig(ys',1,[],xs');
%     m.pv=1*eye(2);  m.pe=diag([0.1 0.001]);
%     xr=[min(y.x); max(y.x)];
%     xhat=pmf(m,y,xr,'Np',6400,'animate','on','delay',0.3)
%     xplot2(y,xhat)

opt=struct('k',0,'Np',100,'animate','off','ind',[],'delay',0,'type','xplot','Xlim',[],'Ylim',[]);
opt=optset(opt,varargin);
if isnan(m.fs)
    error('NL2.PMF: only discrete time models')
end

if ~isnumericscalar(opt.k) | opt.k<0 | opt.k>1
   error('NL2.PMF: k must be 0 or 1')
end
if ~isnumericscalar(opt.Np) | opt.Np<1
   error('NL2.PMF: Np must a positive scalar')
end

if isempty(opt.ind)
   ind=':';
elseif ~isnumeric(opt.ind)
   error('NL2.PMF: value for argument animate must be numeric')
else
   ind=opt.ind;
end
if ~isnumeric(opt.delay) | opt.delay<0
   error('NL2.PMF: delay must be numeric and positive')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
[N,nyz,nuz]=size(z);
if nyz~=ny
   error('NL2.PMF: Number of outputs ny must be the same in model and signal')
end
if nuz~=nu
   error('NL2.PMF: Number of inputs nu must be the same in model and signal')
end
Np=opt.Np;
if nx==1
  xp=linspace(xrange(1,1),xrange(2,1),Np)';
  pv = pdf(m.pv,xp-mean(xp));
elseif nx==2
  N1=round(sqrt(Np));
  Np=N1^2;
  x1=linspace(xrange(1,1),xrange(2,1),N1);
  x2=linspace(xrange(1,2),xrange(2,2),N1);
  [X1,X2]=meshgrid(x1,x2);
  X11=X1(:);
  X22=X2(:);
  xp=[X11 X22];
  pv = pdf(m.pv,xp-ones(Np,1)*mean(xp,1));
  pv=reshape(pv,N1,N1);
else
   error('NL2.PMF: nx must be less than or equal to 2 in the PMF')
end

y=z.y.';
if nu>0;
   u=z.u.';
else
   u=zeros(0,N);
end
t=z.t;
xhat=zeros(N,nx);
Px=zeros(Np,N,nx);
yhat=zeros(N,ny);
Py=zeros(Np,N,ny);
w=ones(Np,1);

for k=1:N;
   % Prediction y(k|k-1)
      try  % vectorized call
         yp=feval(m.h,t(k),xp.',u(:,k).',m.th).';
      catch
         for i=1:Np
            yp(i,:)=feval(m.h,t(k),xp(i,:).',u(:,k),m.th).';
         end
      end
   if opt.k==1
      xhat(k,:)=w'*xp;
      Px(:,k,:)=xp;
      yhat(k,:)=w'*yp;
      Py(:,k,:)=yp;
   end
   % Measurement update of x(k|k)
   w=w.*pdf(m.pe,repmat(y(:,k).',Np,1)-yp);     % Likelihood
%plot([yp,repmat(y(:,k)',Np,1)]),size(yp),pause
   w=w/(eps+sum(w));
%reshape(w,N1,N1),pause
   % Animate
   if strcmp(opt.animate,'on')
      if strcmp(opt.type,'xplot')
         if strcmp(ind,':'), ind=1:nx; end
         for i=1:length(ind)
            subplot(length(ind),1,i)
            %hist(xp(:,ind(i)))
            stem(xp(:,ind(i)),w(:))
            title([char(z.xlabel{ind(i)}),' at time ',num2str(t(k))])
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            if ~isempty(z.x) % True value available
               hold on
               ax=axis;
               plot(z.x(k,ind(i))*[1 1],[ax(3);ax(4)],'r')
               hold off
            end
         end
         drawnow
         pause(opt.delay)
      elseif strcmp(opt.type,'xplot2')
         if strcmp(ind,':'), ind=[1 2]; end
         if length(ind)~=2
            error('NL2.PF: ind must have length 2 in xplot2')
         end
            hp=plot(xp(:,ind(1)),xp(:,ind(2)),'b.');
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            if ~isempty(opt.Ylim), set(gca,'Ylim',opt.Ylim), end
            if ~isempty(z.x) % True value available
               hold on
               plot(z.x(k,ind(1)),z.x(k,ind(2)),'*r')
               hold off
               title(['Time ',num2str(t(k))])
               xlabel(z.xlabel{ind(1)})
               ylabel(z.xlabel{ind(2)})
            end
         drawnow
         pause(opt.delay)
         if k<N
            set(hp,'Xdata',zeros(1,Np),'Ydata',zeros(1,Np))
         end
      elseif strcmp(opt.type,'plot')
         if strcmp(ind,':'), ind=1:ny; end
         for i=1:length(ind)
            subplot(length(ind),1,i)
            %hist(xp(:,ind(i)))
            stem(yp(:,ind(i)),w(:))
            if ~isempty(opt.Xlim), set(gca,'Xlim',opt.Xlim), end
            hold on
            ax=axis;
            plot(z.y(k,ind(i))*[1 1],[ax(3);ax(4)],'r')
            hold off
            title([z.ylabel{ind(i)},' at time ',num2str(t(k))])
         end
         drawnow
         pause(opt.delay)
      end
   end
   if opt.k==0
      % Compute y(k|k)
        try  % vectorized call
           yp=feval(m.h,t(k),xp.',u(:,k).',m.th).';
        catch
           for i=1:Np
              yp(i,:)=feval(m.h,t(k),xp(i,:).',u(:,k),m.th).';
           end
        end
      xhat(k,:)=w'*xp;
      Px(:,k,:)=xp;
      yhat(k,:)=w'*yp;
      Py(:,k,:)=yp;
   end
   % Time update of x(k+1|k)
   try  % vectorized call
      xp=feval(m.f,t(k),xp.',u(:,k),m.th).';
   catch
      for i=1:Np
         xp(i,:)=feval(m.f,t(k),xp(i,:).',u(:,k),m.th).';
      end
   end
   if nx==2
     w=reshape(w,N1,N1);
   end
   w=convn(w,pv,'same');
   w=w(:);
   w=w/(eps+sum(w));
end
%plot(xp(1,:),xp(2,:),'r.')
zhat=sig(yhat,z.t,u.',xhat,Py,Px);
zhat.fs=z.fs;
try
  zhat.xlabel=m.xlabel;
  zhat.ylabel=m.ylabel;
  zhat.ulabel=m.ulabel;
  zhat.name=m.name;
end %PMF
end

function shat=calibrate(s,y,varargin);
%CALIBRATE computes the NL2S parameter estimate in the nl2 (or sensormod) object s from y
%   shat=calibrate(s,y,Property1,Value1,...);
%
%   The states in y.x are assumed known, and the measurements y.y are
%   used to calibrate the parameters in s.th.
%   shat=s except for shat.th and shat.P
%
%   See help nl2s for options on property-value pairs
%
%   Example:
%     %Calibrate a sensormod
%     s=exsensor('toa',3,10);
%     s.th=[0 0 1 0 2 0]; s.pe=0.0001*eye(30);
%     y=simulate(s,1);
%     s0=s; s0.th=s.th+0.5*randn(6,1);
%     shat=calibrate(s0,y);
%     plot(s,s0,shat)

%   Calls nl2.estimate with x0mask=zeros(nx,1)
nx=s.nn(1);
shat=estimate(s,y,'x0mask',zeros(nx,1),varargin{:});
end



function [mhat,res]=estimate(m,z,varargin)
%ESTIMATE estimates/calibrates the parameters in an NL2 system using data
%   [mhat,res]=estimate(m,z,property1,value1,...)
%   The NL2S algorithm is used, which applies the Gauss-Newton algorithm for the
%   parameters and initial states
%
%   m  is an NL2 object
%   z   is a SIG object
%   res is the output struct from the NL2S function
%
%   If multiple datasets z are available, information from each of
%   them is added, so a loop over data is possible using the
%   assignment m=estimate(m,z{i})
%   This works for both x0 and th and any masking
%
%   Property   Value       Description
%   ----------------------------------
%   thmask     ones(nth,1) Masking vector for the parameters
%   x0mask     ones(nx,1)  Masking vector for the initial state
%
%   The optional values are passed to nl2s. See help nl2s for other options.
%
%   See also:
%      calibrate, estimates the parameters, but not x0.
%
%   Examples:
%     %Dynamic system
%     m0=nl2('-th(1)*x^2-th(2)*x','x',[1 0 1 2]);
%     m0.th=[1;1];
%     m0.x0=1;
%     m0.fs=1;
%     z=simulate(m0,0:10);
%     m=m0;           % No model error
%     m.x0=0.8;        % Prior on initial state
%     m.th=[0.9;1.1]; % Prior on parameters
%     mhat=estimate(m,z)
%
%     %Static system (time series)
%     t=(0:0.2:3)';
%     h='th(1)*(1-exp(-th(2)*t))';
%     m=nl2('[]',h,[0 0 1 2]);
%     m.th=[1;1];    % Initial gues
%     m0=m;
%     m0.th=[2;0.5]; % True parameters
%     y=simulate(m0,t)+ndist(0,0.01);  % Noisy data
%     mhat=estimate(m,y)
%     plot(y,simulate(mhat,t),'conf',90)


if nargin<2
   error('NL2 estimate: two input arguments required')
elseif ~isa(m,'nl2')
   error('NL2 estimate: input argument 1 must be an NL2 object')
elseif ~isa(z,'sig')
   if iscell(z)
      for kk=1:length(z)
         if ~isa(z{kk},'sig')
            error('NL2 estimate: input argument 2 must be (an array of) SIG object')
         end
      end
   else
      error('NL2 estimate: input argument 2 must be a SIG object')
   end
end
[mhat,res]=nl2s(m,z,varargin{:});
if ~all(all(m.I==0))
    % Fusion of information in m and mhat
    % Uses that NL2S is ML, not Bayesian, so P from m is not used
    I0=m.I;
    th0=m.th;
    x0=m.x0;
    I1=mhat.I;
    th1=mhat.th;
    x1=mhat.x0;
    I=I1+I0;
    P=pinv(I);
    [E,D]=eig(I);
    d=diag(D);
    ind0=find(d==0);
    ind1=find(d~=0);
    par0=[th0;x0];
    par1=[th1;x1];
    par=par1;
    par(ind1)=P(ind1,ind1)*(I0(ind1,ind1)*par0(ind1) + I1(ind1,ind1)*par1(ind1));
    mhat.th=par(1:length(th0));
    mhat.x0=par(length(th0)+1:end);
    nth=m.nn(5);
    mhat.px0=P(nth+1:end,nth+1:end);
    mhat.P=P(1:nth,1:nth);
    mhat.I=I;
end

mhat.name=[m.name,' (calibrated from data)'];
end


function z=simulate(m,in1,varargin)
%SIMULATE simulates the NL2 system
%   z=simulate(m,z,Property1,Value1,...)
%         with input z.y (SIG object) evaluated in z.t
%   z=simulate(m,T,Property1,Value1,...)
%         otherwise, were T is the simulation time.
%         T may also be a vector for continuous time systems
%
%   1. Continuous time systems
%   The input SIG object u defines the simulation time.
%   If the NL2 system contains no input, use T=[t0 tfinal] instead of u.
%   z is the simulated SIG object, where z.x is the solution to
%     x' = m.f(t, x, u.y; m.th)
%      y = m.h(t, x, u.y; m.th)
%   x(0) = m.x0
%   The initial values x0 and nominal parameters th are the ones defined
%   in the NL2 object.
%   The time instants in z.t are the ones generated by the ODE solver when
%   T=[t0 tfinal], otherwise the output z.y is computed at the time instants
%   specified in T or u.t.
%
%   2. Discrete time systems
%   The system recursions are evaluated in a for loop for t=t0:tfinal
%
%   Property   Value       Description
%   -----------------------------------------------------------------
%   MC         {30}        Number of MC simulations when th uncertain
%
%   Examples:
%     m=exnl('vdp')
%     z=simulate(m,10);
%     plot(z), figure, xplot2(z)


if nargin<2
   error('NL2.simulate: two input arguments required')
end
opt=struct('MC',30);
opt=optset(opt,varargin);

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
if isa(in1,'sig')
   u=in1.y.';
   t=in1.t;
   if ~isnan(in1.fs) & in1.fs~=m.fs
      error('NL2.simulate: m and z must have the same sampling frequency')
   end
   if size(u,1) ~= nu
      error('NL2.simulate: input dimension of m must match output dimension of z')
   end
elseif isvector(in1)
   if isnan(m.fs) % cont time
      if length(in1)==1
          t=[0 in1];
      else
          t=in1;
      end
   elseif length(in1)==1 & nx>0
      t=0:1/m.fs:in1;
   elseif length(in1)==1 & nx==0
      t=in1;  % Static system, one sample
   elseif length(in1)==2
      t=in1(1):1/m.fs:in1(2);
   else
      t=in1;
      Ts=t(2)-t(1);
      if any(abs(diff(t)-Ts)>1e-12)
         error('NL2.simulate: second input as a time vector must have fs consistent with the model')
      end
      if m.fs~=1/Ts;
          error('The sample sampling frequency of m does not correspond to the one in the time vector')
      end
   end
   if nu>0
      u=zeros(nu,length(t));
   else
      u=[]; %XXXzeros(1,length(t));
   end
else
   error('NL2.simulate: second input argument must be a vector or SIG object')
end
N=length(t);
t=t(:);

if isnan(m.fs)  % Continuous time system
   if isempty(m.pv)
      v=zeros(nv,N);
   else
      v=rand(m.pv,N)';
   end
   if nx>0
      if nu>0
         [T,X]=ode45(m.f,[t(1) t(2)],m.x0,[],v(:,1),u(:,1),m.th);
         x=[m.x0 X(end,:).'];
         for k=2:N-1
            [Tnew,Xnew]=ode45(m.f,[t(k) t(k+1)],x(:,end),[],v(:,k),u(:,k),m.th); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T=[T; Tnew(:)];
            X=[X; Xnew];
            x=[x Xnew(end,:).'];
         end
      else
	[T,X]=ode45(m.f,[t(1) t(2)],m.x0,[],[],m.th);
         x=[m.x0 X(end,:).'+v(:,1)];
         for k=2:N-1
            [Tnew,Xnew]=ode45(m.f,[t(k) t(k+1)],x(:,end),[],[],m.th);
            T=[T; Tnew(:)];
            X=[X; Xnew];
            x=[x Xnew(end,:).'+v(:,k)];
         end
      end
      if length(t)==2; % in1 was final time, t is now time vector
         t=T;
         x=X.';
         u=[]; %XXXzeros(nu,length(t));
      end
   else  % No dynamics
      x=[];
   end
else % discrete time
   if nx>0
      x(:,1)=m.x0;
      if isempty(m.pv)
         v=zeros(nv,N);
      else
         v=rand(m.pv,N)';
      end
      if nu>0
         for k=1:N-1
            x(:,k+1)=feval(m.f,t(k),x(:,k),v(:,k),u(:,k),m.th); %%%%%%%%%%%%%%%%
         end
      else
         for k=1:N-1
            x(:,k+1)=feval(m.f,t(k),x(:,k),v(:,k),[],m.th); %%%%%%%%%%%%%%%%%%%%%%%%%
         end
      end
   else
      x=[];
   end
end
y=[];
if nx>0
   if isempty(m.pe)
      e=zeros(ny,N);
   else
      e=rand(m.pe,N)';
   end
   try  % Vectorized definition, which runs faster
      y=feval(m.h,t(:),x,e,u,m.th);
   end
   if size(y,2)~=ny | size(ny,1)~=length(t) % loop over time
      if nu>0
         for k=1:N
             y(:,k)=feval(m.h,t(k),x(:,k),e(:,k),u(:,k),m.th);
         end
      else
         for k=1:N
             y(:,k)=feval(m.h,t(k),x(:,k),e(:,k),[],m.th);
         end
      end
   end
else % static system
      if nu>0
         try  % Vectorized definition, which runs faster
            y=feval(m.h,t(:),[],e,u,m.th);
         catch
            for k=1:length(t);
               y(:,k)=feval(m.h,t(k),[],e(:,k),u(:,k),m.th);
            end
         end
      else
         try
            y=feval(m.h,t(:),[],e,[],m.th);
         catch
            for k=1:length(t);
               y(:,k)=feval(m.h,t(k),[],e(:,k),[],m.th);
            end
         end
      end
      if size(y,1)==ny & size(y,2)==length(t)
         %ok
      elseif size(y,2)==ny & size(y,1)==length(t)
         y=y.';
      else
         error('NL2.SIMULATE: incorrect output from h')
      end
end
% $$$ if isa(m.pe,'pdfclass')
% $$$     y=y+rand(m.pe,size(y,2))';
% $$$ elseif ~isempty(m.pe)
% $$$     y=y+sqrtcov(m.pe)*randn(size(y));
% $$$ end
%if all(u==0); u=[]; end

yMC=[];
xMC=[];
if ~isempty(m.P)
       thMC=rand(ndist(zeros(nth,1),m.P(1:nth,1:nth)),opt.MC)';
       mtmp=m;
       for i=1:opt.MC;
           mtmp.th=m.th+thMC(:,i);
           ztmp=simulate(mtmp,in1,varargin{:},'MC',0);
           yMC(i,:,:)=ztmp.y;
           xMC(i,:,:)=ztmp.x;
       end
end
z=sig(y.',t,u.',x.',yMC,xMC);
z.fs=m.fs;
z.MC=opt.MC;
z.xlabel=m.xlabel;
z.ylabel=m.ylabel;
z.ulabel=m.ulabel;
end

function [mout,zout]=nl22lss(m,z)
%NL22SS returns a linearized state space model
%   Syntax [mout,zout]=nl22ss(m,z)
%
%   The linear state space model is defined by the following Taylor expansion
%      x+  = f(z.t,z.x,z.u) + df(z.t,z.x,z.u)/dx * (x(t)-z.x)
%                           + df(z.t,z.x,z.u)/du * (u(t)-z.u) + v(t)
%     y(t) = h(z.t,z.x,z.u) + dh(z.t,z.x,z.u)/dx * (x(t)-z.x)
%                           + dh(z.t,z.x,z.u)/du * (u(t)-z.u) + e(t)
%   x+ denotes either dx/dt or x(t+1) for cont. and disc. time models, resp.
%   The numeric gradients A=df/dx, B=df/du, C=dh/dx and D=dh/dx are computed
%   around the linearization point specified in the SIG object z (one sample).
%
%   The model is then
%      x+  = ux(t) + A*x(t) + B*u(t) + v(t)
%     y(t) = uy(t) + C*x(t) + D*u(t) + e(t)
%   where the extra inputs are defined as
%      ux(t) = f(z.t,z.x,z.u) - A * z.x - B * z.u
%      uy(t) = h(z.t,z.x,z.u) - C * z.x - D * z.u
%   Using an augmented input vector ua(t)=[u' ux' uy']', the returned model is
%      mout <-> [A, [B I 0], C, [D 0 I]]
%      zout.y = z.y
%      zout.x = z.x
%      zout.u = [z.u; ux; uy];
%
%   Example:
%     mnl2=nl2('-x(1,:).^2-x(1,:)','x',[1 0 1 0])
%     mss=nl22lss(mnl2,sig(0,mnl2.fs,[],2))
%
%   See also: lss.lss2nl2

if nargin<2
   error('NL2 nl22ss: two input arguments required')
elseif ~isa(m,'nl2')
   error('NL2 nl22ss: first input argument must be an NL2 object')
elseif ~isa(z,'sig')
   error('NL2 nl22ss: second input argument must be a SIG object')
elseif size(z,1)>1
   error('NL2 nl22ss: SIG object z has to be one sample, not a signal')
elseif  isempty(z.x)
   error('NL2 nl22ss: SIG object z must contain a state z.x')
end

nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
N=size(z);
h=sqrt(eps);

if size(z.x,2)~=nx
   error('NL2 nl22ss: SIG object state z.x does not match the state dimension in m')
elseif nu>0 & isempty(z.u)
   error('NL2 nl22ss: model m has an input, but the SIG object z.u is empty')
elseif nu~=size(z.u)
   error('NL2 nl22ss: model m input has not the same dimension as the SIG object z.u')
end


A=nl2numgrad(m,z,'dfdx')';
C=nl2numgrad(m,z,'dhdx')';
if nu>0
   B=nl2numgrad(m,z,'dfdu')';
   D=nl2numgrad(m,z,'dhdu')';
   ux = feval(m.f,z.t,z.x.',z.u.',m.th) - A * z.x.' - B * z.u.';
   uy = feval(m.h,z.t,z.x.',z.u.',m.th) - C * z.x.' - D * z.u.';
else
   B=[];
   D=[];
   ux = feval(m.f,z.t,z.x.',z.u.',m.th) - A * z.x.';
   uy = feval(m.h,z.t,z.x.',z.u.',m.th) - C * z.x.';
end


mout=lss(A,[B eye(nx) zeros(nx,ny)],C,[D zeros(ny,nx) eye(ny)],m.fs);

zout=sig(z.y,z.t,[z.u ux.' uy.'],z.x);
zout.fs=z.fs;
end


function der=nl2numgrad(m,z,derstr,mask)
%NL2NUMGRAD returns a numeric gradient approximation
%   Syntax                     Approx  Dim         Application
%   dfdx=nl2numgrad(m,z,'dfdx')   df/dx  (nx,nx,N)    A(t)' in EKF
%   dfdu=nl2numgrad(m,z,'dfdu')   df/du  (nu,nx,N)    B(t)' in EKF
%   dhdx=nl2numgrad(m,z,'dhdx')   dh/dx  (nx,ny,N)    C(t)' in EKF
%   dhdu=nl2numgrad(m,z,'dhdu')   dh/du  (nu,ny,N)    D(t)' in EKF
%   dfdth=nl2numgrad(m,z,'dfdth') df/dth (nth,nx,N)   Estimation
%   dhdth=nl2numgrad(m,z,'dhdth') dh/dth (nth,ny,N)   Estimation
%   dxdth=nl2numgrad(m,z,'dxdth') dx/dth (nth,nx,N)   Estimation
%   dydth=nl2numgrad(m,z,'dydth') dy/dth (nth,ny,N)   Estimation
%   dxdx0=nl2numgrad(m,z,'dxdx0') dx/dx0 (nx,nx,N)    Estimation
%   dydx0=nl2numgrad(m,z,'dydx0') dy/dx0 (nx,ny,N)    Estimation
%   The first four cases correspond to Taylor expansions of the model
%   The last four cases correspond to Taylor expansions of trajectories
%
%   A fourth argument is a binary mask vector for x or th, respectively.

% Prepared for symbolic or user specified gradients, though not implemented
% Adaptive step h in preparation

if nargin<3
   error('NL2 nl2numgrad: at least three input arguments required')
elseif ~isa(m,'nl2')
   error('NL2 nl2numgrad: first input argument must be an NL2 object')
elseif ~isa(z,'sig')
   error('NL2 nl2numgrad: second input argument must be a SIG object')
elseif ~isstr(derstr)
   error('NL2 nl2numgrad: third input argument must be a string')
end

h=sqrt(eps);
nx=m.nn(1);
nv=m.nn(2);
nu=m.nn(3);
ny=m.nn(4);
nth=m.nn(5);
N=size(z);
% Set second input argument to simulate
if nu>0
   in1=sig(z.u,z.t);
   in1.fs=z.fs;
else
   in1=z.t;
end

if strcmpi(derstr,'dfdx') | strcmpi(derstr,'dhdx')
    if nargin<4; mask=ones(nx,1); end
    if sum(mask)~=length(find(mask)),
       error('NL2 nl2numgrad: mask must be binary'),
    end
    n=sum(mask);
    if strcmpi(derstr,'dfdx')
       g=m.f;
       der=zeros(n,nx,N);
    else
       g=m.h;
       der=zeros(n,ny,N);
    end
    Imask=diag(mask);
    indmask=find(mask);
    for j=1:n
      try
         g1=feval(g,z.t,z.x.'-repmat(h*Imask(:,indmask(j)),1,N),z.u.',m.th);
         g2=feval(g,z.t,z.x.'+repmat(h*Imask(:,indmask(j)),1,N),z.u.',m.th);
      catch
         if isempty(z.u); u=zeros(N,1); else, u=z.u; end
         for k=1:N
            g1(:,k)=feval(g,z.t(k),z.x(k,:).'-h*Imask(:,indmask(j)),u(k,:).',m.th);
            g2(:,k)=feval(g,z.t(k),z.x(k,:).'+h*Imask(:,indmask(j)),u(k,:).',m.th);
         end
      end
      der(j,:,:)=(g2-g1)/(2*h);
    end
elseif strcmpi(derstr,'dfdu') | strcmpi(derstr,'dhdu')
    if nargin<4; mask=ones(nu,1); end
    if sum(mask)~=length(find(mask)),
       error('NL2 nl2numgrad: mask must be binary'),
    end
    n=sum(mask);
    if strcmpi(derstr,'dfdu')
       g=m.f;
       der=zeros(n,nx,N);
    else
       g=m.h;
       der=zeros(n,ny,N);
    end
    Imask=diag(mask);
    indmask=find(mask);
    for j=1:n
      try
         g1=feval(g,z.t,z.x.',z.u.'-h*Imask(indmask(j),:),m.th);
         g2=feval(g,z.t,z.x.',z.u.'+h*Imask(indmask(j),:),m.th);
      catch
         for k=1:N
            g1(:,k)=feval(g,z.t(k),z.x(k,:).',z.u(k,:).'-h*Imask(:,indmask(j)),m.th);
            g2(:,k)=feval(g,z.t(k),z.x(k,:).',z.u(k,:).'+h*Imask(:,indmask(j)),m.th);
         end
      end
      der(j,:,:)=(g2-g1)/(2*h);
    end
elseif strcmpi(derstr,'dfdth') | strcmpi(derstr,'dhdth')
    if nargin<4; mask=ones(nth,1); end
    if sum(mask)~=length(find(mask)),
       error('NL2 nl2numgrad: mask must be binary'),
    end
    n=sum(mask);
    if strcmpi(derstr,'dfdth')
       g=m.f;
       der=zeros(n,nx,N);
    else
       g=m.h;
       der=zeros(n,ny,N);
    end
    Imask=diag(mask);
    indmask=find(mask);
    for j=1:n
      try
         g1=feval(g,z.t,z.x.',z.u.',m.th-h*Imask(:,indmask(j)));
         g2=feval(g,z.t,z.x.',z.u.',m.th+h*Imask(:,indmask(j)));
      catch
         if isempty(z.u); u=zeros(N,1); else, u=z.u; end
         for k=1:N
            g1(k,:)=feval(g,z.t(k),z.x(k,:).',u(k,:).',m.th(:)-h*Imask(:,indmask(j)));
            g2(k,:)=feval(g,z.t(k),z.x(k,:).',u(k,:).',m.th(:)+h*Imask(:,indmask(j)));
         end
      end
      der(j,:,:)=(g2-g1)/(2*h);
    end
elseif strcmpi(derstr,'dxdx0') | strcmpi(derstr,'dydx0')
    if nargin<4; mask=ones(nx,1); end
    if sum(mask)~=length(find(mask)),
        error('NL2 nl2numgrad: mask must be binary'),
    end
    n=sum(mask);
    if strcmpi(derstr,'dxdx0')
       der=zeros(n,nx,N);
    else
       der=zeros(n,ny,N);
    end
    Imask=diag(mask);
    indmask=find(mask);
    mtmp=m;
    mtmp.pe=[];  % Denoisify simulation model
    for j=1:n
       mtmp.x0=m.x0-h*Imask(:,indmask(j));
       z1=simulate(mtmp,in1);
       mtmp.x0=m.x0+h*Imask(:,indmask(j));
       z2=simulate(mtmp,in1);
       if strcmpi(derstr,'dxdx0')
          der(j,:,:)=(z2.x-z1.x).'/(2*h);
       else
          der(j,:,:)=(z2.y-z1.y).'/(2*h);
       end
    end
elseif strcmpi(derstr,'dxdth') | strcmpi(derstr,'dydth')
    if nargin<4; mask=ones(nth,1); end
    if sum(mask)~=length(find(mask)),
        error('NL2 nl2numgrad: mask must be binary'),
    end
    n=sum(mask);
    if strcmpi(derstr,'dxdth')
       der=zeros(n,nx,N);
    else
       der=zeros(n,ny,N);
    end
    Imask=diag(mask);
    indmask=find(mask);
    mtmp=m;
    mtmp.pe=[];  % Denoisify simulation model
    for j=1:n
       mtmp.th=m.th-h*Imask(:,indmask(j));
       z1=simulate(mtmp,in1);
       mtmp.th=m.th+h*Imask(:,indmask(j));
       z2=simulate(mtmp,in1);
       if strcmpi(derstr,'dxdth')
          der(j,:,:)=(z2.x-z1.x).'/(2*h);
       else
          der(j,:,:)=(z2.y-z1.y).'/(2*h);
       end
    end
else
     error('NL2 nl2numgrad: unknown option for der')
end
end


end % methods
end
