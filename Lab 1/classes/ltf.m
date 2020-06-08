classdef ltf  < ltimod
%LTF Transfer function implementation of Linear Time Invariant (LTI) systems
%   Causal MIMO deterministic systems supported.
%
%   ltf.a, ltf.b         Polynomials in transfer function H(q) = b(q)/a(q)
%                                B(q)
%                        y(t)=  ----- u(t-nk) +
%                                A(q)
%   nn                   Polynomial orders [na,nb,nk,nu,ny]
%   ltf.fs               Sampling frequency
%
%   help ltf.ltf      gives help for the constructor
%   methods ltf       lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
   a, b, nn
   fs, name, ulabel, ylabel, tlabel, marker, markerlabel, desc
   MC, sysMC
end

methods
% =================================
% -----------Constructor-----------
% =================================

function m=ltf(varargin)
%LTF  Construct a Transfer Function
%            B(q)
%     y(t)= ------ u(t-nk)
%            A(q)
%   or
%     a(1)y(t)+a(2)y(t-1)+...+a(na)y(t-na+1)=
%          b(1)u(t)+b(2)u(t-1)+...+b(nb)u(t-nb+1)
%
%   To define a LTF object, use
%      sys=ltf(b,a,fs)     % Discrete time LTF, fs in Hertz
%      sys=ltf(b,a)        % Continuous time by default, with fs=NaN convention
%   To define a structure for estimation or random model generation, use
%      sys=ltf([na,nb,nk])          % SISO
%      sys=ltf([na,nb,nk,nu,ny])    % MIMO
%   where na, nb are the orders of a and b, respectively, and nk is the delay.
%
%   fs=NaN for continuous time systems
%   The entered numerator b and denominator a polynomials can be arbitrary
%   vectors of arbitrary lengths. However, it is recommended to fill up with
%   zeros to equal length. LTF will otherwise fill up with zeros FROM THE LEFT.
%   This corresponds to polynomials in descending powers of s and z, resp.
%
%   Inside the object, b and a are normalized using the following conventions:
%   1. Coefficient a(1) preceding y(t) is always one.
%   2. b(1) is non-zero
%   3. If b and a are specified of different length, the shorter one is
%      extended with zeros from the right.
%      Use b,a of equal size to avoid problems!
%   4. nk is the relative degree, so nk>=0 for causal systems,
%      and nk<0 for non-causal systems
%   That is, the stored difference equation (or differential equation
%   for continuous time) is
%     y(t)+a(1)y(t-1)+...+a(na)y(t-na)=
%         b(1)u(t-nk)+b(2)u(t-nk-1)+...+b(nb)u(t-nk-nb+1)
%
%   For MIMO systems, the SISO transfer function from ui to yj is
%   given by b(j,:,i) and a.
%   Note that a and nk are the same for all inputs and outputs, in order
%   to be consistent with (unstructured) state space models.
%
%
%   Special constructions:
%   s=ltf('unit')        % Unit function
%   s=ltf('s')           % Laplace operator
%   q=ltf('q',fs)        % Time shift operator
%   z=ltf('z',fs)        % z transform operator (same as q)
%   qinv=ltf('delay')    % Inverse time shift
%
%   Examples:
%     Gc=ltf([1 0 0],[1 1 1])
%     Gc=ltf(1,[1 1 1])        % Note convention!
%     Gd=ltf(1,[1 1 1],2)


a=[]; b=[];
fs=NaN;
MC=0;
sysMC=[];
name='';  xlabel=[]; ulabel=[]; ylabel=[]; desc=[];


if nargin<1
    error('LTF contructor expects at least one input argument: see help ltf')
elseif isa(varargin{1},'ltf')  % LTF
    s=varargin{1};
    if nargin>1
       MC=varargin{2};
    else
       MC=30;
    end
    s2=s;
    % Unpack
    b=s2.b; a=s2.a; fs=s2.fs; nn=s2.nn;
    name=s2.name; ulabel=s2.ulabel; ylabel=s2.ylabel; desc=s2.desc;
elseif isa(varargin{1},'lss')     % Conversion from SS object
    s=varargin{1};
    s2=lss2ltf(s);
    % Compute MC samples
    if ~isempty(s.sysMC)
        MC=length(s.sysMC);
        srand=rand(s,MC);             %MC samples of s...
        for k=1:MC
            sysMC{k}=ltf(srand{k});     %are converted to ss
        end
    end
    % Unpack
    b=s2.b; a=s2.a; fs=s2.fs; nn=s2.nn;
    name=s2.name; ulabel=s2.ulabel; ylabel=s2.ylabel; desc=s2.desc;
elseif isa(varargin{1},'ltf')  % LSS
        s2=lss2ltf(s);
        % Compute MC samples
        if ~isempty(s.sysMC)
            MC=length(s.sysMC);
            srand=rand(s,MC);             %MC samples of s...
            for k=1:MC
                sysMC{k}=ltf(srand{k});     %are converted to ss
            end
        else
            MC=0;
        end
    % Unpack
    b=s2.b; a=s2.a; fs=s2.fs; nn=s2.nn;
    name=s2.name; ulabel=s2.ulabel; ylabel=s2.ylabel; desc=s2.desc;
elseif isa(varargin{1},'arx')    % Conversion from ARX object
    s2=arx2ltf(s);
    % Compute MC samples
    if ~isempty(s.P) & ~all(s.P==0);
        if size(s.P,1)~=size(s.P,2)
           MC=size(s.P,2);
        elseif s.MC>0
           MC=s.MC;
        else
           MC=30;
        end
        srand=rand(s,MC);             %MC samples of s...
        for k=1:MC
           sysMC{k}=ltf(srand{k});     %are converted to ss
        end
    else
        MC=0;
    end
    % Unpack
    b=s2.b; a=s2.a; fs=s2.fs; nn=s2.nn;
    name=s2.name; ulabel=s2.ulabel; ylabel=s2.ylabel; desc=s2.desc;
elseif isstr(varargin{1}) % Special LTF
    if strcmpi(varargin{1},'unit')
        b=1; a=1; nn=[1 1 0 1 1]; fs=NaN;
    elseif strcmpi(varargin{1},'delay') | strcmpi(varargin{1},'qinv')
        b=1; a=1; nn=[0 1 1 1 1]; fs=1;
    elseif strcmpi(varargin{1},'s')
        b=1; a=1; nn=[1 1 -1 1 1]; fs=NaN;
    elseif strcmpi(varargin{1},'q')
        b=1; a=1; nn=[1 1 -1 1 1]; fs=1;
    elseif strcmpi(varargin{1},'z')
        b=1; a=1; nn=[1 1 -1 1 1]; fs=1;
    else
        error('LTF: Unknown string option')
    end
    if nargin>1
        error(sprintf('''%s'' may not be followed by any arguments.', ...
              varargin{1}))
    end
elseif nargin==1 & isvector(varargin{1}) % empty structure
    nn=varargin{1};
    if ~isnumeric(nn) | size(nn,1)>1
       error('LTF constructor: nn must be a numeric row vector in a call ltf(nn)')
    end
    na=nn(1);
    if length(nn)>1; nb=nn(2); else, nb=na; end
    if length(nn)>2; nk=nn(3); else, nk=0; end
    if length(nn)>3; nu=nn(4); else, nu=1; end
    if length(nn)>4; ny=nn(5); else, ny=1; end
    if length(nn)>5
      error('LTF constructor: Too long [na,nb,nk,nu,ny] vector.')
    end
    nn=[na nb nk nu ny];
    a=[]; b=[];
elseif nargin==3 | nargin==2;
    b=varargin{1};
    a=varargin{2};
    if ~isnumeric(a) | ~isnumeric(b)
      error('LTF constructor: a and b must be numeric in ltf(b,a).')
    end
    if isempty(a) | isempty(b)
      error('LTF constructor: a and b cannot be empty in ltf(b,a).')
    end
    if nargin>2,  fs=varargin{3}; else, fs=NaN; end
    if fs==0 | isnan(fs)
        fs=NaN;
    end
    if ~isnumericscalar(fs) | fs<=0
        error('LTF constructor: fs must be positive or NaN')
    end
    ny=size(b,1);
    nb=size(b,2);
    nu=size(b,3);
    na=size(a,2);
    if nb<na   % Convention: Fill up with initial zeros in b
        b=[zeros(ny,na-nb,nu) b];
        nb=na;
    end
    if nb>na   % Convention: Fill up with initial zeros in a
        a=[zeros(1,nb-na) a];
        na=nb;
    end
    % Remove initial and trailing zeros
    inda=find(abs(a)>1e-14); % Remove zeros with rounding errors
    a=a(inda(1):inda(end));
    indb=zeros(1,nb); % All non-zero positions
    for i=1:nu
        for j=1:ny
            indb=indb | (abs(b(j,:,i))>1e-14);
        end
    end
    indb=find(indb);  % Index with non-zero coefficients in at least one channel
    if isempty(indb)
        %error('All parameters in b are zero')
        b=zeros(ny,1,nu);
        indb=inda;
    else
        b=b(:,indb(1):indb(end),:);
    end
    % normalize
    b=b/a(1);
    a=a/a(1);
    % Determine delay and orders
    nk=indb(1)-inda(1);
    na=length(a)-1;
    nb=size(b,2);
    nn=[na nb nk nu ny];
else
    disp('Usage: sys=ltf(b,a)                      % Contnuous time')
    disp('       sys=ltf(b,a,fs)                   % Discrete time')
    disp('       sys=ltf([na,nb,nk])               % SISO structure')
    disp('       sys=ltf([na,nb,nk,ny,nu])         % MIMO structure')
    error('LTF constructor: Incorrect number of inputs')
end

m.a=a; m.b=b;
m.fs=fs; m.nn=nn;
m.name=name; m.ulabel=ulabel; m.ylabel=ylabel; m.desc=desc;

end %constructor


% =================================
% -----Non-static functions--------
% =================================


function out=fieldread(arg)
out=eval(arg);
end

function fieldwrite(arg1,arg2)
if isnumeric(arg2) | isa(arg2,'cell')
   if strcmp(arg1,'b')
      if isa(arg2,'cell')
         error('b cannot be assiged a cell')
      end
      if isequal(size(arg2),size(b))
         b=arg2;
         %eval([arg1,'=',mat2str(arg2),';'])
      else
         error('Cannot change size of b in LTF objects')
      end
   elseif strcmp(arg1,'a')
      if isa(arg2,'cell')
         error('b cannot be assiged a cell')
      end
      if isequal(size(arg2),size(a))
         a=arg2;
      else
         error('Cannot change size of a in LTF objects')
      end
   elseif strcmp(arg1,'sysMC')
      if isempty(arg2)
         sysMC=[];
         MC=0;
      elseif ~isa(arg2,'cell')
         error('sysMC must be a cell array')
      elseif ~isa(arg2{1},'ltf')
         error('sysMC must be a cell array of LTF objects')
      elseif isequal(size(arg2{1}.nn),nn)  % Check sysMC argument
         sysMC=arg2;
      else
         sysMC=arg2;
         %error('Argument of sysMC does not match the LTF dimensions')
      end
   else
      error(['There is no visible field called ' arg1 '.'])
   end
else
   error('Field value must be numeric, cell or PDFCLASS.')
end
end %function

function sys=arrayread(jj,ii)
%ARRAYREAD used to pick out sub-systems by indexing
%   sys=arrayread(jj,ii)
%   jj is the row index/indices, corresponding to the outputs
%   ii is the column index/indices, corresponding to the inputs
%   Example: s(2,3) gives the SISO system from input 3 to output 2

if nargin<2
   error('LTF.ARRAYREAD: both input indeces in G(j,i) are required')
end
nk=nn(3);
nu=nn(4);
ny=nn(5);
if isstr(ii)
    if strcmp(ii,':')
        ii=1:nu;
    else
       error('LTF.ARRAYREAD: Inputs ii must be a vector or :')
    end
end
if isstr(jj)
    if strcmp(jj,':')
        jj=1:ny;
    else
       error('LTF.ARRAYREAD: Outputs jj must be a vector or :')
    end
end
if any(ii<=0) | any(ii~=round(ii))
     error('LTF.ARRAYREAD: Input (column) index must be positive integers')
end
if any(jj<=0) | any(jj~=round(jj))
     error('LTF.ARRAYREAD: Output (column) index must be positive integers')
end
if any(ii>nu)
     error('LTF.ARRAYREAD: Input (column) index larger than the number of inputs')
end
if any(jj>ny)
     error('Output (row) index larger than the number of outputs')
end
nkfill=nk+size(b,2)-size(a,2);
btmp=cat(2,b(jj,:,ii),zeros(length(jj),-nkfill,length(ii)));
atmp=[a zeros(1,nkfill)];
sys=ltf(btmp,atmp,fs);
stmp2=[];
for k=1:length(sysMC)
   stmp1=sysMC{k};
   stmp2{k}=stmp1(jj,ii);
end
sysMC=stmp2;

end %function


function sys=uplus(sys)
%UPLUS unitary plus
%   sys=uplus(sys) or sys=+sys
sys=sys;

end %function


function sys=uminus(sys)
%UMINUS  unitary minus
%   sys=uminus(sys) or sys=-sys
sys=sys;
sys.b=-sys.b;

end %function


function sys=minus(s1,s2)
%MINUS parallel connection with difference at output
%   sys=minus(s1,s2) or sys=s1-s2
%   Requirements: nu1=nu2, ny1=ny2

sys=s1+(-s2);

end %function


function sys=power(s,n)
%POWER is repeated multiplication
%   sys=power(s,n)  or  s.^n
%   calls s^n
sys=s^n;

end %function


function sys=mpower(s,n)
%MPOWER is repeated multiplication
%   sys=mpower(s,n)  or  s^n

if ~isnumeric(n) || ~isscalar(n) || round(n)~=n
    error('Power must be an integer')
end
if n==0;
    sys=ltf('unit');
elseif n<0
    sys=mpower(inv(s),-n);
else
    sys=s;
    for k=2:n
        sys=sys*s;
    end
end

end %function


function sys=inv(s)
%INV computes the right inverse of a system such that s*sys=I
%   sys=inv(s)
%   Requirements: ny=nu=1 for non-causal LTF, only SISO inverse implemented currently
%   Otherwise, ny<=nu, where ss/inv is invoked

nu=s.nn(4);
ny=s.nn(5);
if nu==1 & ny==1
    sys=ltf(s.a,s.b,s.fs);
    sys.nn(3)=-s.nn(3);
elseif ny<=nu
    sys=ltf(inv(lss(s)));
else
    error('ny<=nu is required for right inverse')
end
end %function


%------- Functions of two systems-------------




function out=eq(s1,s2)
%EQ equality
%   out=eq(s1,s2) or s1==s2

[s1,s2]=ltfcheck(s1,s2);
[b1,a1]=unpack(s1);
[b2,a2]=unpack(s2);
% NaN==NaN is false
%s1.nn,s2.nn,b1,b2,a1,a2
out=false;
if all(s1.nn==s2.nn)
    if size(b1)==size(b2) & size(a1)==size(a2)
        if all(b1-b2<1e-10) & all(a1-a2<1e-10) %& s1.fs==s2.fs
            out=true;
        end
    else
        disp('Warning: cancellation might occur, try minreal first')
    end
end
end %function


function [b,d,h]=isequal(s1,s2,conf);
%ISEQUAL performs a hypothesis test with confidence level conf
%   [b,d,h]=isequal(s1,s2,conf);
%   b is true or false
%   d is the chi2 test statistic or [] if the structures of s1 and s2 differ
%   h is the threshold implied by conf
%   conf is the confidence level (default 0.99)
%   See also: eq (non-stochastic version)

if nargin<3, conf=0.99; end
if any(s1.nn~=s2.nn)
    b=false; d=[]; h=[];
else
    m1=arx(s1);
    m2=arx(s2);
    [b,d,h]=isequal(m1,m2,conf);
end
end %function


function sys=mrdivide(s1,s2)
%MRDIVIDE computes the right inverse of a system such that s*sys=I
%   sys=s1*inv(s2) is used
[s1,s2]=ltfcheck(s1,s2);
sys=s1*inv(s2);
sys.sysMC=mceval2('mrdivide',s1,s2);
end %function


function sys=mldivide(s1,s2)
%MLDIVIDE computes the left inverse of a system such that s*sys=I
%   sys=inv(s1)*s2 is used
[s1,s2]=ltfcheck(s1,s2);
sys=inv(s1)*s2;
sys.sysMC=mceval2('mldivide',s1,s2);
end %function


function sys=plus(s1,s2)
%PLUS Add two models together in a parallel connection
%   sys=plus(s1,s2)  or sys=s1+s2
%   Requirements: nu1=nu2, ny1=ny2
%   Special treatment if s1 or s2 is a scaling matrix.
%   If s1 or s2 is a matrix, element (i,j) gives the direct term
%   from input i to output j

[s1,s2]=ltfcheck(s1,s2);
if isa(s1,'ltf') | isa(s2,'ltf')
    if isa(s1,'ltf') & isa(s2,'ltf')
        fs=fscheck(s1,s2);
    end
    if isnumeric(s1)
        b1(:,1,:)=s1; a1=1;  ny1=size(b1,1); nu1=size(b1,3);
    elseif isa(s1,'ltf')
        [b1,a1]=unpack(s1);
        fs=s1.fs;
        nu1=s1.nn(4); ny1=s1.nn(5);
    else
        error('s1 must be either a LTF object or a matrix')
    end
    if isnumeric(s2)
        b2(:,1,:)=s2; a2=1;   ny2=size(b2,1); nu2=size(b2,3);
    elseif isa(s2,'ltf')
        [b2,a2]=unpack(s2);
        fs=s2.fs;
        nu2=s2.nn(4); ny2=s2.nn(5);
    else
        error('s2 must be either a LTF object or a matrix')
    end
    if ny1~=ny2 | nu1~=nu2
       disp(['nu1=',num2str(nu1),', ny1=',num2str(ny1),', nu2=',num2str(nu2),', ny2=',num2str(ny2)])
       error('Incorrect dimensions in s1+s2. nu1=nu2, ny1=ny2 required.')
    end
    for i=1:nu1;
        for j=1:ny1;
            b(j,:,i)=conv(b1(j,:,i),a2)+conv(b2(j,:,i),a1);
        end
    end
    a=conv(a1,a2);
    sys=ltf(b,a,fs);
    sys.sysMC=mceval2('plus',s1,s2);
end
end %function


function sys=mtimes(s1,s2)
%MTIMES (*) multiplies two models together in a series connection
%   sys=mtimes(s1,s2)  or sys=s1*s2
%   Implements the series connection y=s1*s2*u
%   Requirements: ny2=nu1
%   Special treatment if s1 or s2 is a (scaling) matrix.
%   If s1 or s2 is a matrix, element (i,j) gives the scaling
%   from input i to output j.

[s1,s2]=ltfcheck(s1,s2);
if isa(s1,'ltf') | isa(s2,'ltf')
    if isa(s1,'ltf') & isa(s2,'ltf')
        fs=fscheck(s1,s2);
    end
    if isnumeric(s1)
        b1(:,1,:)=s1; a1=1;   ny1=size(b1,1); nu1=size(b1,3);
    elseif isa(s1,'ltf')
        [b1,a1]=unpack(s1);
        fs=s1.fs;
        nu1=s1.nn(4); ny1=s1.nn(5);
    else
        error('s1 must be either a LTF object or a matrix')
    end
    if isnumeric(s2)
        b2(:,1,:)=s2; a2=1;   ny2=size(b2,1); nu2=size(b2,3);
    elseif isa(s2,'ltf')
        [b2,a2]=unpack(s2);
        fs=s2.fs;
        nu2=s2.nn(4); ny2=s2.nn(5);
    else
        error('s2 must be either a LTF object or a matrix')
    end
    if ny2~=nu1
       disp(['ny2=nu1 is required in mtimes, but ny2=',num2str(ny2),' and nu1=',num2str(nu1)])
       error('Incorrect dimensions in s1*s2.')
    end
    b=zeros(ny1,size(b1,2)+size(b2,2)-1,nu2);
    for i=1:nu2;
        for j=1:ny1;
            for k=1:nu1; %=ny2, matrix multiplication
                b(j,:,i)=b(j,:,i)+conv(b2(k,:,i),b1(j,:,k));
            end
        end
    end
    a=conv(a1,a2);
    sys=ltf(b,a,fs);
    sys.sysMC=mceval2('mtimes',s1,s2);
end
end %function


function sys=times(s1,s2)
%TIMES (.*) multiplies two models together elementwise
%   sys=times(s1,s2)  or sys=s1.*s2
%   Implements the connection y=(s1.*s2)*u
%   Requirements: nu1=nu2 and ny1=ny2
%   Special treatment if s1 or s2 is a (scaling) matrix.
%   If s1 or s2 is a matrix, element (i,j) gives the scaling
%   from input i to output j.

[s1,s2]=ltfcheck(s1,s2);
if isa(s1,'ltf') | isa(s2,'ltf')
    if isa(s1,'ltf') & isa(s2,'ltf')
        fs=fscheck(s1,s2);
    end
    if isnumeric(s1)
        b1(:,1,:)=s1; a1=1;   ny1=size(b1,1); nu1=size(b1,3);
    elseif isa(s1,'ltf')
        [b1,a1]=unpack(s1);
        fs=s1.fs;
        nu1=s1.nn(4); ny1=s1.nn(5);
    else
        error('s1 must be either a LTF object or a matrix')
    end
    if isnumeric(s2)
        b2(:,1,:)=s2; a2=1;   ny2=size(b2,1); nu2=size(b2,3);
    elseif isa(s2,'ltf')
        [b2,a2]=unpack(s2);
        fs=s2.fs;
        nu2=s2.nn(4); ny2=s2.nn(5);
    else
        error('s2 must be either a LTF object or a matrix')
    end
    if nu1~=nu2
       error(['LTF: Incorrect dimensions in s1.*s2, nu1=nu2 is required but nu1=',num2str(nu1),' and nu2=',num2str(nu2)])
    end
    if ny1~=ny2
       error(['LTF: Incorrect dimensions in s1.*s2, ny1=ny2 is required but ny1=',num2str(ny1),' and ny2=',num2str(ny2)])
    end
    b=zeros(ny1,size(b1,2)+size(b2,2)-1,nu1);
    for i=1:nu1;
        for j=1:ny1;
            b(j,:,i)=b(j,:,i)+conv(b2(j,:,i),b1(j,:,i));
        end
    end
    a=conv(a1,a2);
    sys=ltf(b,a,fs);
end
end %function


function sys=feedback(s1,s2)
%FEEDBACK computes the feedback connection of two systems
%   sys=feedback(s1,s2)
%   System s1 in forward loop and s2 in feedback loop.
%   Requirements: nu1=ny2, nu2=ny1

[s1,s2]=ltfcheck(s1,s2);
if isa(s1,'ltf') | isa(s2,'ltf')
    if isa(s1,'ltf') & isa(s2,'ltf')
        fs=fscheck(s1,s2);
    end
    if isnumeric(s1)
        b1(:,1,:)=s1; a1=1;   ny1=size(b1,1); nu1=size(b1,3);
    elseif isa(s1,'ltf')
        [b1,a1]=unpack(s1);
        fs=s1.fs;
        nu1=s1.nn(4); ny1=s1.nn(5);
    else
        error('s1 must be either a LTF object or a matrix')
    end
    if isnumeric(s2)
        b2(:,1,:)=s2; a2=1;   ny2=size(b2,1); nu2=size(b2,3);
    elseif isa(s2,'ltf')
        [b2,a2]=unpack(s2);
        fs=s2.fs;
        nu2=s2.nn(4); ny2=s2.nn(5);
    else
        error('s2 must be either a LTF object or a matrix')
    end
    if  ny2>1
        error('Currently, only scalar feedback path (s2*s1 scalar) is implemented for LTF')
    end
    if ny1~=nu2
       disp(['ny1=nu2 and ny2=nu1 are required, '])
       disp(['but nu1=',num2str(nu1),' nu2=',num2str(nu2),' ny1=',num2str(ny1),' and ny2=',num2str(ny2)])
       error('Incorrect dimensions in feedback(s1,s2).')
    end
    g0=s2*s1;
    sys=inv(eye(g0)+g0)*s1;
    sys=inherit(sys,s1,'Closed loop system');
    sys.sysMC=mceval2('feedback',s1,s2);
end
end %function


function sys=diag(varargin)
%DIAG appends independent models into one using append recursively
%   sys=diag(s1,s2,s3,...)  % Observe, no brackets
%   sys=diag(s)  with one input argument extract the diagonal of the
%   transfer function matrix s

if nargin==1;
   s=ltf(varargin{1});
   n=min(s.nn(4:5));
   sys=s(1,1);
   for k=2:n
      sys=diag(sys,s(k,k));
   end
   sys.sysMC=mceval1('diag',varargin{1});
elseif nargin>1;
   sys=append(ltf(varargin{1}),varargin{2});
   sys.sysMC=mceval2('append',ltf(varargin{1}),varargin{2});
end
if nargin>2  % Recursion
   sys=diag(sys,varargin{3:end});
   sys.sysMC=mceval2('diag',sys,varargin{3:end});
end
end %function


function sys=append(s1,s2)
%APPEND packs two independent models into one.
%   sys = append(sys1,sys2)
%   Requirements: none

[s1,s2]=ltfcheck(s1,s2);
if isa(s1,'ltf') | isa(s2,'ltf')
    if isa(s1,'ltf') & isa(s2,'ltf')
        fs=fscheck(s1,s2);
    end
    if isnumeric(s1)
        b1(:,1,:)=s1'; a1=1;   ny1=size(b1,1); nu1=size(b1,3);
    elseif isa(s1,'ltf')
        [b1,a1]=unpack(s1);
        fs=s1.fs;
        nu1=s1.nn(4); ny1=s1.nn(5);
    else
        error('s1 must be either a LTF object or a matrix')
    end
    if isnumeric(s2)
        b2(:,1,:)=s2'; a2=1;   ny2=size(b2,1); nu2=size(b2,3);
    elseif isa(s2,'ltf')
        [b2,a2]=unpack(s2);
        fs=s2.fs;
        nu2=s2.nn(4); ny2=s2.nn(5);
    else
        error('s2 must be either a LTF object or a matrix')
    end
    b=zeros(ny2,size(b1,2)+size(b2,2)-1,nu1);
    for i=1:nu1;
        for j=1:ny1;
                b1new(j,:,i)=conv(b1(j,:,i),a2);
        end
    end
    for i=1:nu2;
        for j=1:ny2;
                b2new(j,:,i)=conv(b2(j,:,i),a1);
        end
    end
    b=zeros(ny1+ny2,size(b1new,2),nu1+nu2);
    for k=1:size(b2new,2);
        b(:,k,:)=[squeeze(b1new(:,k,:)) zeros(ny1,nu2);zeros(ny2,nu1) squeeze(b2new(:,k,:))];
    end
    a=conv(a1,a2);
    sys=ltf(b,a,fs);
    sys.sysMC=mceval2('append',s1,s2);
end
end %function


% =================================
% ----------Concatenations---------
% =================================



function sys=ctranspose(s)
%CTRANSPOSE reverse the inputs with outputs for MIMO systems
%   sys=ctranspose(s)  or sys=s'
%   ctranspose is the same as transpose
sys=transpose(s);
sys.sysMC=mceval1('transpose',s);

end %function


function sys=transpose(s)
%TRANSPOSE reverse the inputs with outputs for MIMO systems
%   sys=transpose(s)  or sys=s.'
%
%   For transfer functions this corresponds symbolically to
%   b -> b'  along the first and third dimensions

[b,a]=unpack(s);
[na,nb,nk,nu,ny]=size(s);
for i=1:size(b,3);
    for j=1:size(b,1);
        bb(i,:,j)=b(j,:,i);
    end
end
sys=ltf(bb,a,s.fs);
sys.sysMC=mceval1('transpose',s);
end %function


function sys=horzcat(varargin)
%HORZCAT performs horizontal concatenation
%   sys=horzcat(s1,s2,...)
%   This can be used to create a MISO system from SISO systems
%   Reguirement: All systems must have the same number of outputs ny1=ny2=...
%   nk is the smallest of all delays
%   The transfer functions are put on a common denominator

if nargin<1;
    error('No system input argument')
elseif nargin<2;
    sys=varargin{1};
else
    s1=varargin{1};
    s2=varargin{2};
    if isa(s1,'ltf') | isa(s2,'ltf')
        if isa(s1,'ltf') & isa(s2,'ltf')
            fs=fscheck(s1,s2);
        end
        if isnumericscalar(s1)
            b1=s1; a1=1;  ny1=1; nu1=1;
        elseif isa(s1,'ltf')
            [b1,a1]=unpack(s1);
            fs=s1.fs;
            nu1=s1.nn(4); ny1=s1.nn(5);
        else
            error('Input must be either a LTF object or a matrix')
        end
        if isnumericscalar(s2)
            b2=s2; a2=1; ny2=1; nu2=1;
        elseif isa(s2,'ltf')
            [b2,a2]=unpack(s2);
            fs=s2.fs;
            nu2=s2.nn(4); ny2=s2.nn(5);
        else
            error('Input must be either a LTF object or a matrix')
        end
        if ny1~=ny2
            error(['Incorrect dimensions in horzcat. ny1=ny2 required, but ny1=',...
                    num2str(ny1),' and ny2=',num2str(ny2)])
        end
        for i=1:nu1;
            for j=1:ny1;
                b(j,:,i)=conv(b1(j,:,i),a2);
            end
        end
        for i=1:nu2;
            for j=1:ny1;
                b(j,:,nu1+i)=conv(b2(j,:,i),a1);
            end
        end
        a=conv(a1,a2);
        sys=ltf(b,a,fs);
    end
end
if nargin>2 % Recursion
   sys=horzcat(sys,varargin{3:end});
end
sys.sysMC=mceval1('horzcat',varargin{:});
end %function


function sys=vertcat(varargin)
%VERTCAT performs vertical concatenation
%   sys=vertcat(s1,s2,...)
%   This can be used to create a SIMO system from SISO systems
%   Reguirement: All systems must have the same number of inputs nu1=nu2=...
%   nk is the smallest of all delays
%   The transfer functions are put on a common denominator


if nargin<1;
    error('No system input argument')
elseif nargin<2;
    sys=varargin{1};
else
    s1=varargin{1};
    s2=varargin{2};
    if isa(s1,'ltf') | isa(s2,'ltf')
        if isa(s1,'ltf') & isa(s2,'ltf')
            fs=fscheck(s1,s2);
        end
        if isnumericscalar(s1)
            b1=s1; a1=1;  ny1=1; nu1=1;
        elseif isa(s1,'ltf')
            [b1,a1]=unpack(s1);
            fs=s1.fs;
            nu1=s1.nn(4); ny1=s1.nn(5);
        else
            error('Input must be either a LTF object or a matrix')
        end
        if isnumericscalar(s2)
            b2=s2; a2=1; ny2=1; nu2=1;
        elseif isa(s2,'ltf')
            [b2,a2]=unpack(s2);
            fs=s2.fs;
            nu2=s2.nn(4); ny2=s2.nn(5);
        else
            error('Input must be either a LTF object or a matrix')
        end
        if nu1~=nu2
            error(['Incorrect dimensions in horzcat. nu1=nu2 required, but nu1=',...
                    num2str(nu1),' and nu2=',num2str(nu2)])
        end
        for i=1:nu1;
            for j=1:ny1;
                b(j,:,i)=conv(b1(j,:,i),a2);
            end
        end
        for i=1:nu2;
            for j=1:ny2;
                b(ny1+j,:,i)=conv(b2(j,:,i),a1);
            end
        end
        a=conv(a1,a2);
        sys=ltf(b,a,fs);
    end
end
if nargin>2 % Recursion
   sys=vertcat(sys,varargin{3:end});
end
sys.sysMC=mceval1('vertcat',varargin{:});
end %function



% =================================
% ----------LTI conversions--------
% =================================



function sys=ltf2lss(s,form)
%LTF2LSS converts a LTF to SS object
%   sys=ltf2lss(s,form)
%
%   s is the LTF object, and sys is a SS object
%   form='observer'   or 'o' for observability form
%   form='controller' or 'c' for controllability form
%
%   Unlike the ltf2lss function, this method works for MIMO.
%   The observability form works straightforwardly for SIMO systems, and
%   the controllability form works for MISO systems.
%   For MIMO systems, there is no simple standard form.
%   Here, a simple append is used, which leads to a non-minimal realization.
%   For instance, for the observability form, one SIMO realization is found for
%   each output, and these are then appended.
%   minreal can be used to decrease the model order, but then the structure is
%   is lost.
%
%   Examples:
%
%   See also: lss.lss2ltf, ltf2lss, lss2ltf

if s.nn(3)<0
    error('LTF.LTF2LSS: Only causal LTF can be converted to SS')
end
if isempty(s.a) & isempty(s.b)
   sys=lss([s.nn(1)+s.nn(3) s.nn(4) 0 s.nn(5)]);
   return
end
[b,a,fs]=unpack(s);
[na,nb,nk,nu,ny]=size(s);
na=length(a);
if size(b,2)>na
   error('LTF.LTF2LSS: Only causal LTF can be converted to SS')
end
if nargin<2,
    if nu==1
        form='controller';
    elseif ny==1
        form='observer';
    else
        form='observer';
    end
end
switch lower(form)
  case {'c', 'controller'}
    % Controller canonical form
    C=[]; D=[];
    A=kron(eye(nu),[-a(2:end); eye(na-2,na-1)]);
    B=kron(eye(nu),eye(na-1,1));
    for i=1:nu
        for j=1:ny
            Ctmp(j,:) = b(j,2:na,i) - b(j,1,i) * a(2:end);
            Dtmp(j,1) = b(j,1,i);
        end
        C=[C Ctmp];
        D=[D Dtmp];
    end
  case {'o', 'observer'}
    % Observer canonical form
    A = kron(eye(ny),[-a(2:end)' eye(na-1,na-2)]);
    C = kron(eye(ny),eye(1,na-1));
    B=[]; D=[];
    for j=1:ny
        for i=1:nu
            Btmp(:,i) = b(j,2:na,i)' - b(j,1,i) * a(2:end)';
            Dtmp(1,i) = b(j,1,i);
        end
        B=[B;Btmp];
        D=[D;Dtmp];
    end
  otherwise
    error(['Illegal form: ' form ])
end
sys=lss(A,B,C,D,fs);
sys=inherit(sys,s,'LTF -> SS');
sys.sysMC=mceval1('ltf2lss',s);
end %function


function sys=minreal(s)
%MINREAL computes the minimal realization
%   sys=minreal(s)
%   The function is based on the output from zpk
%   First, common zeros and poles are cancelled by systematically
%   looping through z and p
%   Then, inputs and outputs that are decoupled are removed by searching
%   the gain matrix k for all zero rows and columns

[z,p,k]=zpk(s);
nu=s.nn(4);
ny=s.nn(5);
if s.nn(1)==0 | s.nn(2)==0
    sys=s;
    return % Nothing to do
end
removepind=[];
for pind=1:length(p);
    ok=1;
    for i=1:nu;
        for j=1:ny
            ind=find(abs(p(pind)-z(j,:,i))<1e-5);
            if isempty(ind);
                ok=0;
            else
                removezind(j,i)=ind(1);
            end
            %common=intersect(common,z(j,:,i));
        end
    end
    if ok
        removepind=[removepind pind];
        for i=1:nu;
            for j=1:ny
                z(j,removezind(j,i),i)=0;
            end
        end
    end
end
p(removepind)=zeros(size(removepind));
a=poly(p);
if isempty(z)
    b=1;
else
    ny=size(z,1);
    nu=size(z,3);
    for i=1:nu;
        for j=1:ny
            ind1=find( isinf(z(j,:,i)));
            ind2=find( ~isinf(z(j,:,i)));
            b1=zeros(size(ind1));
            b2=poly(z(j,ind2,i));
            b(j,:,i)=[b1 b2];
        end
    end
end
nk=s.nn(3);
ii=[];
jj=[];
for j=1:size(k,1);  % Loop through outputs to find decoupled outputs
    if any(k(j,:)~=0)
       jj=[jj j];
    end
end
for i=1:size(k,2);  % Loop through inputs to find decoupled inputs
    if any(k(:,i)~=0)
       ii=[ii i];
    end
end
%btmp=cat(2,b(jj,:,ii),zeros(length(jj),-nk,length(ii)));
%atmp=[a zeros(1,nk)];
%sys=k(jj,ii).*ltf(btmp,atmp,s.fs);
sys=k(jj,ii).*ltf(b(jj,:,ii),a,s.fs);
sys=inherit(sys,s,'Minimal realization');
sys.sysMC=mceval1('minreal',s);
end %function


function sys=c2d(s,fs,method)
%C2D converts continuous time LTF to discrete time LTF
%   sys=c2d(s,fs,method)
%
%   fs {1} is the sampling frequency
%   method is a string describing the assumption on intersample behaviour
%   method is a string describing the assumption on intersample behaviour
%   {'ZOH'}    zero order hold, piece-wise constant input assumed
%   'FOH'      first order hold, piecewise linear input assumed
%   'bilinear' s=2/T (z-1)/(z+1)
%
%   Example:
%     Gc=getfilter(6,0.2,'fs',NaN);
%     Gd=c2d(Gc,1)
%     zpplot(Gc,Gd)
%
%   See also:

if nargin<3; method='zoh'; end
if nargin<2; fs=1; end
if s.fs>0
    error('LTF object already in discrete time')
end
T=1/fs;

method=lower(method);
if strcmp(method(1:3),'bil');
    [b,a]=unpack(s);
    ny=size(b,1);
    nu=size(b,3);
    na=length(a);
    pm{1}=[1]; pp{1}=[1];
    for i=1:na-1
        pm{i+1}=conv(pm{i},[1 -1]);
        pp{i+1}=conv(pp{i},[1 1]);
    end
    for i=1:na;
        W(i,:)=(2/T)^(na-i)*conv(pm{na+1-i},pp{i});
    end
    a=a*W;
    for j=1:ny
        for i=1:nu
            b(j,:,i)=b(j,:,i)*W;
        end
    end
    b=b/a(1); a=a/a(1);
    sys=ltf(b,a,fs);
elseif strcmp(method(1:3),'zoh');
    sys=ltf(c2d(lss(s),fs,'zoh'));
elseif strcmp(method(1:3),'foh');
    sys=ltf(c2d(lss(s),fs,'foh'));
else
    error(['LTF.C2D: Method ',method,' unknown'])
end
sys=inherit(sys,s,'Discretized model');
sys.sysMC=mceval1('c2d',s);
end %function


function sys=d2c(s,method)
%D2C converts discrete time LTF to continuous time LTF
%   sys=d2c(s,method)
%
%   method is a string describing the assumption on intersample behaviour
%   {'ZOH'}    zero order hold, piece-wise constant input assumed
%   'FOH'      first order hold, piecewise linear input assumed
%   'bilinear' s=2/T (z-1)/(z+1),  z=-(s+2/T)/(s-2/T)

if nargin<2; method='zoh'; end
if isnan(s.fs)
    error('LTF object already in continuous time')
end

T=1/s.fs;

method=lower(method);
if strcmp(method(1:3),'bil');
    [b,a]=unpack(s);
    na=length(a);
    pm{1}=[1]; pp{1}=[1];
    for i=1:na-1
        pm{i+1}=conv(pm{i},[-1 -2/T]);
        pp{i+1}=conv(pp{i},[1 -2/T]);
    end
    for i=1:na;
        W(i,:)=conv(pm{na+1-i},pp{i});
    end
    a=a*W;
    b=b*W;
    b=b/a(1); a=a/a(1);
    sys=ltf(b,a);
elseif strcmp(method(1:3),'zoh');
    sys=ltf(d2c(lss(s),'zoh'));
elseif strcmp(method(1:3),'foh');
    sys=ltf(d2c(lss(s),'foh'));
else
    error(['Method ',method,' unknown'])
end
sys=inherit(sys,s,'Continuous model computed from discrete model');
sys.sysMC=mceval1('d2c',s);
end %function


function m=ltf2arx(s,varargin)
%LTF2ARX computes the ARX correspondence of a LTF object
m=arx(s.b,s.a);
m.pe=0;
if ~isempty(s.sysMC)
   MC=length(s.sysMC);
   for k=1:MC;
      mtmp=arx(s.sysMC{k}.b,s.sysMC{k}.a);
      P(:,k)=mtmp.th;
   end
else
   P=[];
end
m=arx(s.nn,m.th,P);
end %function


function Hf=ltf2ft(s,varargin)
%LTF2FT computes frequency domain response H(f) of a LTF object
%   Hf=ltf2ft(ltf,Property1,Value1,...)
%
%   N number of linearly spaced grid points on frequency axis
%   between 0 and fs/2
%   Fields in ft
%   f   frequency grid points
%   H   frequency response H(f,j,i) for output j and input i
%   HMC Monte Carlo samples
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   N         {1024}      Number of frequency grid points
%   fmax      {'auto'}    Maximum frequency
%   f         {}          Frequency grid (overrides N and fmax)
%
%   Examples:
%     G=rand(ltf([6 6 0 2 2]));
%     Gf=ltf2ft(G); % Explicit call
%     Gf=ft(G);    % Implicit call
%     plot(Gf)
%
%   See also:

opt=struct('N',1024,'f',[],'fmax','auto','MC',0);
opt=optset(opt,varargin);
fs=s.fs;
[b,a]=unpack(s);
[na,nb,nk,nu,ny]=size(s);

if isempty(opt.f);
    if isnan(fs);  % Continuous time ltf
        if isnumericscalar(opt.fmax) & opt.fmax>0
            fmax=opt.fmax;
        else
            T=timeconstant(s.a,fs);
            if isempty(T); T=10; end %ad-hoc
            fmax=8/T;
            % ad-hoc rule: Maximum interesting frequency is 8 times larger than
            % dominating pole and zero
        end
        f=linspace(0,1,opt.N)*fmax;
    else
        f=linspace(0,0.5,opt.N)*fs;
    end
else
    f=opt.f;
    f=f(:).';
end
N=length(f);
n=length(a);
if fs>0;
    W=exp(-i*(0:n-1)'*2*pi*f/fs);
else
    W= repmat((i*2*pi*f),n,1).^repmat((n-1:-1:0)',1,N);
end
for j=1:ny
    for l=1:nu
        H(:,j,l) = ((b(j,:,l)*W)./(a*W)).';
    end
end
HMC=[];
MC=length(s.sysMC);
for k=1:MC
    Htmp=freq(s.sysMC{k});
    HMC(k,:,:,:) = Htmp.H;
end
Hf=ft(H,f,HMC);
Hf=inherit(Hf,s,'LTF -> FREQ');
end %function


function [z,p,k,zMC,pMC,kMC]=zpk(s)
%ZPK computes the zeros, poles and gain
%   [z,p,k,zMC,pMC,kMC]=zpk(s)
%   z is a (ny x nb x nu) matrix with zeros
%     for strictly proper transfer functions, z is filled up with Inf
%     according to the system theory definition G(Inf)=0 for proper systems
%   p is a row vector with poles
%   k is a (ny x nu) matrix with gains
%   zMC,pMC,kMC are the corresponding Monte Carlo arrays
%
%   Examples:
%     G=rand(ltf([2 2 0 2 2]))
%     [z,p,k]=zpk(G)
%   See also: lss.zpk

[b,a]=unpack(s);
nk=s.nn(3);
nu=s.nn(4);
ny=s.nn(5);
z=[];
for i=1:nu;
    for j=1:ny
        ind=find(b(j,:,i)~=0);
        if isempty(ind)
            nktmp=length(b(j,:,i))-1;
        else
            nktmp=ind(1)-1;  % Initial zeros in b(j,:,i)
        end
        %nktmp=size(z,2);  % better way to fill up: NO
        broots=[roots(b(j,:,i)).' Inf*ones(1,nktmp)];
        if ~isempty(broots)
           z(j,:,i)=broots;
        else
           z(j,:,i)=zeros(1,size(z,2),1);
        end
        k(j,i)=b(j,nktmp+1,i);
    end
end
p=[roots(a).']; % zeros(1,nk)];
if ~isempty(s.sysMC)
   for i=1:length(s.sysMC)
      [zMC(i,:,:,:),pMC(i,:),kMC(i)]=zpk(s.sysMC{i});
   end
else
   pMC=[];
   zMC=[];
   kMC=[];
end
end %function


% =================================
% ----------LTI generators---------
% =================================



function sys=eye(s)
%EYE generates a unitary transfer function matrix of the same size as s.
%    Example: sys=eye(ltf([1,1,0,2,2]));

[na,nb,nk,nu,ny]=size(s);
if nu~=ny
    error('EYE works only for square systems')
else
    b(1:nu,1,1:nu)=eye(nu);
    sys=ltf(b,1,s.fs);
end
end %function


function sys=zeros(s)
%ZEROS generates a zero transfer function matrix of the same size as s.
%    Example: sys=zeros(ltf([1,1,0,2,2]));

[na,nb,nk,nu,ny]=size(s);
sys=ltf(zeros(ny,1,nu),1,s.fs);
end %function


function sys=ones(s)
%ONES generates a zero transfer function matrix of the same size as s.
%    Example: sys=ones(ltf([1,1,0,2,2]));

[na,nb,nk,nu,ny]=size(s);
sys=ltf(ones(ny,1,nu),1,s.fs);
end %function


% Stochastic systems


function Gu=uncertain(G,c,X,MC)
%UNCERTAIN introduces uncertainty in coefficients c according to dist X
%   Gu=uncertain(G,c,X,MC)
%
%   c   is a string containing the coefficients to randomize,
%       e.g. c='a(2), b(1)  b(3)' (the separators are arbitrary).
%   X   is a distribution of the same length (length(X))
%   MC  is the number of MC samples (default G.MC)
%   Gu  is the same as G, except for that the field sysMC has changed
%       and the mean values of the coefficients in c are set to E(X)
%
%   Example:
%   G=ltf([1 2 3],[1 1.2 1],1);
%   Gu=uncertain(G,'a(2)',udist(1.1,1.3),100)

if nargin<4; MC=G.MC; end
if ~isnumeric(MC) | length(MC)>1 | MC<0 | MC~=round(MC)
    error('LTF.UNCERTAIN: MC must be a positive integer')
elseif MC==0
   Gu=G;
   return
end

if isa(X,'pdfclass')
   ind=find(c=='(');
   ind2=find(c==')');
   N=length(ind);
   if N==length(X)
      Gu=G;                % Copy
      Gu.sysMC=rand(G,MC); % assure that GuMC has correct dim
      Gu.MC=MC;
      x=rand(X,MC);        % random coefficients
      mu=E(X);
      for k=1:N
         cstr{k}=c(ind(k)-1:ind2(k));  % coefficient string
         eval(['Gu.',cstr{k},'=mu(k);'])
         for m=1:MC          % randomize the coefficients in sysMC
            eval(['Gu.sysMC{m}.',cstr{k},'=x(m,k);'])
         end
      end
   else
      error('LTF.UNCERTAIN: when setting a number of coefficients to a PDFCLASS, the length of X must match')
   end
else
   error('LTF.UNCERTAIN: X must be a PDFCLASS object')
end
end %function


function sys=rand(s,MC,varargin)
%RAND generates a random LTF
%   sys=rand(s,MC,options)
%
%   sys is a cell array with MC samples if MC>1.
%   Default MC=1, in case sys is a LTF object (no cell array)
%   sys is empty if s has no uncertainty
%
%   1. If s is an empty model with only structure specified, then
%      the random models are taken from a general prior distribution.
%      The function randpoly is used to generate poles and zeros,
%      and the options are forwarded to randpoly
%   2. If s is an uncertain model, then MC samples from this are generated.
%      The MC samples of the model are taken from the field sysMC, and
%      these may in turn have been obtained from direct definition of
%      (a) uncertainty in the constructor, (b) estimation uncertainty,
%      or (c) Monte Carlo simulation uncertainty in estimation
%   3. If G is a LTF model without uncertainty, then MC copies are returned
%
%   Example:
%      rand(ltf([2,2,1,2,1]) generates a random model or order 2 with
%                        2 inputs and 1 output
%      rand(sys) generates a random model of the same size as sys
%                if sys is uncertain.
%
%   See also: lss.rand, randpoly

[na,nb,nk,nu,ny]=size(s);
fs=s.fs;
if nargin<2; MC=1; end
if isempty(s.sysMC) & isempty(s.a)  % Models from the prior
    nktmp=nk-na+nb-1;
    for k=1:MC;
        [a,fs]=randpoly(na,'fs',fs,varargin{:});
        a=[a zeros(1,nktmp)];
        if nu==0; b=1; end
        for i=1:nu;
            for j=1:ny;
                b(j,:,i)=[randpoly(nb-1,'fs',fs,varargin{:}) zeros(1,-nktmp)];
            end
        end
        sys{k} = ltf(b,a,fs);
    end
elseif length(s.sysMC)>0  % MC samples
    if MC<=length(s.sysMC)
        ind=1:MC; % Systematic sampling
    else
        ind=ceil(length(s.sysMC)*rand(1,MC));  % Sampling with replacement
    end
    for k=1:MC;
        sys{k} = s.sysMC{ind(k)};
    end
else  % MC replicas
    for k=1:MC
       sys{k}=s;
    end
end
if MC==1 & ~isempty(sys), sys=sys{1}; end
end %function


function sys=E(s)
%E returns the expected value of the system by removing MC data
%   sys=E(s)
if isempty(s.a)  % Empty structure
    sys=ltf(s.nn);
else
    sys=ltf(s.b,s.a,s.fs);
end
end %function


function s=fix(su)
%FIX removes uncertainty in LTF coefficients
%   s=fix(su)
%
%   This might be useful to speed up computations by avoiding
%   the Monte Carlo loop over uncertain parameters
%
%   See also: lss.fix, ltf.uncertain, ltf.estimate

s=su;
s.sysMC=[];
end %function


function Gout=estimate(Gin,z,varargin)
%ESTIMATE estimates a state space model from data in a SIG object
%   Gout=estimate(Gin,z,Property1,Value1,...)
%
%   A LTF model with structure as specified in Gin is estimated from
%   the signal z using a two-step least squares (LS) algorithm.
%   1. A high order FIR model is estimated using LS as a first step.
%   2. In the second step, the high order FIR model is simulated without noise.
%   3. Then the low order ARX model is estimated using LS.
%   4. Finally, the ARX model is converted into a LTF with the same
%      b and a polynomials.
%
%   For model uncertainty representation, it is assumed that
%      a. the true system is contained in the high order FIR model
%      b. and that the estimate can be considered Gaussian distributed.
%   The latter is  true for Gaussian noise and asymptotically otherwise.
%   Monte Carlo simulations are used to convert the Gaussian high order
%   estimate to a sample based representation of the non-Gaussian distribution
%   of the low order ARX estimate.
%
%
%   Property  Value/Default   Description
%   MC        {30}            Number of Monte Carlo simulations
%   nfir      min([10*na,50]) FIR order in the first step
%
%
%   Example:
%     G=c2d(rand(ltf(2)));
%     y=filter(G,getsignal('prbs',100))+0.1*randn(100,1);
%     Ghat=estimate(ltf(2),y);
%     plot(G,Ghat)
%
%   See also: arx.estimate, nl.estimate


[na,nb,nk,nu,ny]=size(Gin);
[N,nysig,nusig,nx]=size(z);
if nu~=nusig
    error('LTF.estimate: The system and signal must have the same number of inputs')
end
if ny~=nysig
    error('LTF.estimate: The system and signal must have the same number of outputs')
end
if ~isempty(z.yMC)
   MC=length(z.yMC);
else
   MC=30;
end
opt=struct('MC',MC,'nfir',min([10*na,50]));
opt=optset(opt,varargin);

% Step 1: high-order FIR estimation
Gtmp1=estimate(arx([0 opt.nfir nk nu ny]),z);


% Step 2: noise-free simulation
Gtmp1.pe=0;
u=sig(z.u,z.fs);
ztmp=simulate(Gtmp1,u);

% Step 3
Gtmp2=estimate(arx([na nb nk nu ny]),ztmp);

% Step 4
Gtmp2.MC=opt.MC;
Gout=ltf(Gtmp2);
end %function


% =================================
% -------Data generators-----------
% =================================


function y=impulse(varargin)
%IMPULSE generates the pulse/impulse response of a system
%   y=impulse(G,T)
%   The analytical solution in lss.impulse is used
%   by calling y=impulse(lss(G),T)
%
%   Argument    Description
%   --------------------------------------------------------
%   G           LTF object
%   T           Simulation length (number of samples N or time span),
%               or a time vector t
%               Default it is  estimated from dominating pole
%   y           Output SIG object
%
%   Examples:
%     Gd=getfilter(4,0.3,'fs',2);
%     yd=impulse(Gd);
%     Gc=getfilter(4,0.3);
%     yc=impulse(Gc);
%     subplot(2,1,1), staircase(yd)
%     subplot(2,1,2), plot(yc)
%
%   See also: lss.impulse, ltf.step, ltf.simulate

G=varargin{1};
if G.nn(3)<0
   error('LTF.IMPULSE only accepts causal systems')
end
if nargout>0
   y=impulse(lss(G),varargin{2:end});
else
   impulse(lss(G),varargin{2:end});
end
return
end %function


function y=step(varargin)
%STEP generates the step response of a system
%   y=step(G,T)
%   The analytical solution in lss.step is used
%   by calling y=step(lss(G),T)
%
%   Argument    Description
%   --------------------------------------------------------
%   G           LTF object
%   T           Simulation length (number of samples N or time),
%               or a time vector t
%               Default it is  estimated from dominating pole
%   y           Output SIG object
%
%   Examples:
%     Gd=getfilter(4,0.3,'fs',2);
%     yd=step(Gd);
%     Gc=getfilter(4,0.3);
%     yc=step(Gc);
%     subplot(2,1,1), staircase(yd)
%     subplot(2,1,2), plot(yc)
%
%   See also: lss.step, ltf.impulse, ltf.simulate

G=varargin{1};
if G.nn(3)<0
   error('LTF.STEP only accepts causal systems')
end
if nargout>0
   y=step(lss(G),varargin{2:end});
else
   step(lss(G),varargin{2:end});
end
end %function


function y=filter(G,u,varargin)
%FILTER filters a signal u with the LTF filter object G to get y=Gu
%   y=filter(G,u,varargin)
%
%   The low-level filter function is used internally.
%   Basically, the following is performed
%     ytmp=filter(G.b,G.a,u.y);
%     y=sig(ytmp,u.fs);
%   Note that the sampling frequency of the input SIG object has precedence
%   to the one specified in G.
%
%   The filter method, unlike the filter function, works for
%   MIMO LTF and SIG objects.
%
%   Monte Carlo data are propagated according to the following precedence rules:
%   1. If the signal contains MC data (u.MC>0), then y gets the same number of
%      Monte Carlo samples, each one corresponding to a filtering to one input
%      realization.
%   2. Otherwise, if the LTF object is uncertain, then the input u is filtered
%      through G.MC Monte Carlo realizations of G.
%
%   Examples:
%     G=getfilter(4,0.3);
%     u=getsignal('prbs');
%     y=filter(G,u);
%     Monte Carlo:
%     u=getsignal('square',256,128);
%     fs=u.fs;
%     u=u.y;
%     uMC=repmat(u',MC,1)+0.1*randn(MC,length(u));
%     u=sig(u,fs,[],[],uMC);
%     y=filter(G,u);
%     staircase(y,'conf',90)
%     MIMO:
%     Gstruc=ltf([2 2 0 2 2]);
%     Gstruc.fs=1;
%     G2=rand(Gstruc)
%     u=getsignal('square',64);
%     u2=[u sig(zeros(size(u),1))]
%     y2=filter(G2,u2);
%     staircase(y2)
%
%   See also: filter, ltf.filtfilt, ltf.ncfilt

fs=u.fs;
uout=u.y;
[N,nu]=size(uout);
[na,nb,nk,nuG,ny]=size(G);
[b,a,Gfs]=unpack(G);
if isempty(a)
   error('LTF.FILTER: The LTF object is empty')
end
if nu~=nuG;
   error('LTF.FILTER: The number of inputs in G (nu) and number of outputs in u (ny) must be equal')
end
if isinf(fs)
    error('LTF.FILTER works only for discrete time models, use simulate otherwise')
end
if isinf(Gfs)
    error('LTF.FILTER works only for discrete time signals, use simulate otherwise')
end

for j=1:ny
   yout(:,j)=filter(G.b(j,:,1),G.a,uout(:,1));
   for i=2:nu
      yout(:,j)=yout(:,j)+filter(G.b(j,:,i),G.a,uout(:,i));
   end
end
yMC=[];
uMC=[];
if ~isempty(u.yMC)  % MC simulations of input
   uMC=u.yMC;
   Gfix=E(G); % Remove uncertainty if any
   for k=1:u.MC;
       utmp=shiftdim(uMC(k,:,:),1);
       ytmp=filter(Gfix,sig(utmp,fs));
       yMC(k,:,:)=ytmp.y;
   end
elseif G.MC>0  % MC simulations of LTF realizations
   for k=1:G.MC;
       ytmp=filter(G.sysMC{k},u);
       yMC(k,:,:)=ytmp.y;
   end
end

y=sig(yout,fs,uout,[],yMC);
y=inherit(y,u,'Filtered signal');
end %function


function y=filtfilt(G,u,varargin)
%FILTFILT filters a signal u forward and backward in time
%   y=filtfilt(G,u,varargin)
%   The LTF filter object G gives with the filtfilt procedure
%       Y(z)=G(z)G(1/z)U(z)=|G(z)|^2U(z)
%   The low-level filtfilt function is used internally.
%   The filtfilt method, unlike the filtfilt function, works for
%   MIMO LTF and SIG objects.
%
%   Basically, the following is performed
%     ytmp=filtfilt(G.b,G.a,u.y);
%     y=sig(ytmp,u.fs);
%   Note that the sampling frequency of the input SIG object has precedence
%   to the one specified in G.
%
%   Monte Carlo data are propagated according to the following precedence rules:
%   1. If the signal contains MC data (u.MC>0), then y gets the same number of
%      Monte Carlo samples, each one corresponding to a filtering to one input
%      realization.
%   2. Otherwise, if the LTF object is uncertain, then the input u is filtered
%      through G.MC Monte Carlo realizations of G.
%
%   Examples:
%     G=getfilter(4,0.3);
%     u=getsignal('prbs');
%     y=filtfilt(G,u);
%     Monte Carlo:
%     u=getsignal('square',256,128);
%     fs=u.fs;
%     u=u.y;
%     uMC=repmat(u',MC,1)+0.1*randn(MC,length(u));
%     u=sig(u,fs,[],[],uMC);
%     y=filtfilt(G,u);
%     staircase(y,'conf',90)
%     MIMO:
%     Gstruc=ltf([2 2 0 2 2]);
%     Gstruc.fs=1;
%     G2=rand(Gstruc)
%     u=getsignal('square',64);
%     u2=[u sig(zeros(size(u),1))]
%     y2=filtfilt(G2,u2);
%     staircase(y2)
%
%   See also: filtfilt, ltf.filter, ltf.ncfilter

fs=u.fs;
uout=u.y;
[N,nu]=size(uout);
[na,nb,nk,nuG,ny]=size(G);
[b,a,Gfs]=unpack(G);

if nu~=nuG;
   error('The number of inputs in G (nu) and number of outputs in u (ny) must be equal')
end
if isinf(fs)
    error('LTF filtfilt method works only for discrete time models, use simulate otherwise')
end
if isinf(Gfs)
    error('LTF filtfilt method works only for discrete time signals, use simulate otherwise')
end

for j=1:ny
   yout(:,j)=filtfilt(G.b(j,:,1),G.a,uout(:,1));
   for i=2:nu
      yout(:,j)=yout(:,j)+filtfilt(G.b(j,:,i),G.a,uout(:,i));
   end
end
yMC=[];
uMC=[];
if u.MC>0 & ~isempty(u.yMC) % MC simulations of input
   uMC=u.yMC;
   Gfix=E(G); % Remove uncertainty if any
   for k=1:u.MC;
       utmp=shiftdim(uMC(k,:,:),1);
       ytmp=filtfilt(Gfix,sig(utmp,fs));
       yMC(k,:,:)=ytmp.y;
   end
elseif G.MC>0  % MC simulations of LTF realizations
   for k=1:G.MC;
       ytmp=filtfilt(G.sysMC{k},sig(u,fs));
       yMC(k,:,:)=ytmp.y;
   end
end

y=sig(yout,fs,uout,[],yMC);
y=inherit(y,u,'Filtered signal (filtfilt)');
end %function


function y=ncfilter(G,u,varargin)
%NCFILTER filters a signal u forward and backward in time
%   y=ncfilter(G,u,varargin)
%   The LTF filter object G gives with the ncfilter procedure
%       Y(z)=Gc(z)Gn(1/z)U(z)
%   where Gc is the causal factor of G with all poles and zeros inside
%   the unit circle and Gn is the non-causal part, which is implemented
%   backwards in time.
%   Basically, the following is performed
%     ytmp=ncfilter(G.b,G.a,u.y);
%     y=sig(ytmp,u.fs,u.y);
%   Note that the sampling frequency of the input SIG object has precedence
%   to the one specified in G.
%
%   The low-level ncfilter function is used internally.
%   The ncfilter method, unlike the ncfilter function, works for
%   MIMO LTF and SIG objects.
%
%   Monte Carlo data are propagated according to the following precedence rules:
%   1. If the signal contains MC data (u.MC>0), then y gets the same number of
%      Monte Carlo samples, each one corresponding to a filtering to one input
%      realization.
%   2. Otherwise, if the LTF object is uncertain, then the input u is filtered
%      through G.MC Monte Carlo realizations of G.
%
%   Examples:
%     G=getfilter(4,0.3);
%     u=getsignal('prbs');
%     y=ncfilter(G,u);
%     Monte Carlo:
%     u=getsignal('square',256,128);
%     fs=u.fs;
%     u=u.y;
%     uMC=repmat(u',MC,1)+0.1*randn(MC,length(u));
%     u=sig(u,fs,[],[],uMC);
%     y=ncfilter(G,u);
%     staircase(y,'conf',90)
%     MIMO:
%     Gstruc=ltf([2 2 0 2 2]);
%     Gstruc.fs=1;
%     G2=rand(Gstruc)
%     u=getsignal('square',64);
%     u2=[u sig(zeros(size(u),1))]
%     y2=ncfilter(G2,u2);
%     staircase(y2)
%
%   See also: ncfilter, ltf.filter, ltf.filtfilt

fs=u.fs;
uout=u.y;
[N,nu]=size(uout);
[na,nb,nk,nuG,ny]=size(G);
[b,a,Gfs]=unpack(G);

if nu~=nuG;
   error('The number of inputs in G (nu) and number of outputs in u (ny) must be equal')
end
if isinf(fs)
    error('LTF ncfilter method works only for discrete time models, use simulate otherwise')
end
if isinf(Gfs)
    error('LTF ncfilter method works only for discrete time signals, use simulate otherwise')
end

for j=1:ny
   yout(:,j)=ncfilter(G.b(j,:,1),G.a,uout(:,1));
   for i=2:nu
      yout(:,j)=yout(:,j)+ncfilter(G.b(j,:,i),G.a,uout(:,i));
   end
end
yMC=[];
uMC=[];
if u.MC>0  & ~isempty(u.yMC) % MC simulations of input
   uMC=u.yMC;
   Gfix=E(G); % Remove uncertainty if any
   for k=1:u.MC;
       utmp=shiftdim(uMC(k,:,:),1);
       ytmp=ncfilter(Gfix,sig(utmp,fs));
       yMC(k,:,:)=ytmp.y;
   end
elseif G.MC>0  % MC simulations of LTF realizations
   for k=1:G.MC;
       ytmp=ncfilter(G.sysMC{k},sig(u,fs));
       yMC(k,:,:)=ytmp.y;
   end
end

y=sig(yout,fs,uout,[],yMC);
y=inherit(y,u,'Filtered signal (ncfilter)');
end %function


function [z,xf]=simulate(G,u,varargin)
%SIMULATE simulates a signal from a LTF model
%   [z,xf]=simulate(G,u,Property1,Value1,...)
%
%   G is the LTF object
%   u is signal object with the same
%   sampling interval as the model (for instance, use u=sig(uvec,fs)).
%   If there are inconsistent sampling intervals, the one in the SIG object
%   has precedence.
%   xf is the final state, and can be forwarded as argument 'xi'
%   for simulation over different segments.
%
%   1. For discrete time systems, the low-level filter functions is applied
%      to each input-output MIMO channel.
%   2. For continuous time systems, the simulation is based on fast sampling.
%      The time constant Tc of the system is computed using timeconstant. This
%      is based on the dominating stable pole. The continuous time system is
%      then resampled 200 times faster using c2d(G,200/Tc). Alternatively,
%      you can set the property Ts at your choice.
%      For signals with discontinuities, these are first located.
%      These can be either steps or impulses, see the SIG object contructor
%      how these work. The input is then segmented between the boundaries
%      defined by the discontinuities, and a separate sampling and simulation
%      is done in each segment, where the filter state is saved and used
%      in the next segment.
%
%
%   For non-uniformly sampled inputs, only continuous models apply.
%
%   Property   Value/{Default}   Description
%   --------------------------------------------------------------------------
%   MC                           Number of Monte Carlo simulations
%                                Default value inherited from model or signal.
%   xi         zeros(nx,1)       Initial state as a nx=max([na nb-1]) vector
%   Ts         {timeconstant(s)/200} Sampling interval for fast sampling
%
%   Example:
%     Gc=exlti('tf3c');
%     fs=5;
%     Gd=c2d(Gc,fs);
%     ud=[getsignal('zeros',5);getsignal('ones',45)];
%     ud.fs=fs;
%     yd=simulate(Gd,ud);
%     umat=[0 0 2 2]';
%     t=[0 1 1 10];
%     u=sig(umat,t);
%     y=simulate(Gc,u);
%     subplot(2,1,1), staircase(yd)
%     subplot(2,1,2), plot(y)
%
%   See also: lss.simulate

if nargin<2;
    error('LTF.SIMULATE: An input must be specified')
end
if ~isa(u,'sig'),
    error('LTF.SIMULATE: The input u must be a SIG object')
end

if ~isempty(u.yMC)
   MCdef=size(u.yMC,1);
elseif ~isempty(G.sysMC)
   MCdef=length(G.sysMC);
else
   MCdef=30;
end

[na,nb,nk,nu,ny]=size(G);
[b,a,fs]=unpack(G);
if a(1)==0
   error('LTF.SIMULATE: non-causal system (a(1)~=0) cannot be simulated')
end
xi=zeros(ny,max([na nb-1]),nu);

opt=struct('MC',MCdef,'xi',xi,'Ts',0);
opt=optset(opt,varargin);

if fs>0
    %-------------------------------
    %---discrete time model-------
    %-------------------------------
    if isnan(u.fs)
        error('LTF.SIMULATE: Signal must be discrete time when model is discrete time')
    end
    umat=u.y;
    [N,nudata]=size(umat);
    if N<nudata; error('LTF.SIMULATE: Input matrix u should have more rows than columns'), end
    if nudata~=nu
       error('LTF.SIMULATE: Number of inputs in u not compatible with model')
    end
    if abs(u.fs-fs)>1e-12
       fs=u.fs;  % Precedence rule
    end
    t=(0:N-1)'/fs;
    for i=1:nu
        for j=1:ny
            [ymat(:,j),xf(j,:,i)]=filter(b(j,:,i),a,umat(:,i),xi(j,:,i));
        end
    end
    % Ready with nominal realization, now create MC realizations
    if opt.MC>0  % user can force it to zero by property value spec
        MC=opt.MC;
        if ~isempty(u.yMC)
           % Monte Carlo simulations from data realizations
           utmp=rand(u,MC);
           for k=1:MC
               ztmp=simulate(G,utmp{k},'MC',0,'Ts',opt.Ts);
               ymatMC(k,:,:)=ztmp.y;
           end
        elseif MC>0 & ~isempty(G.sysMC)
           % Monte Carlo simulations from model realizations
           Gtmp=rand(G,MC);
           for k=1:MC
               ztmp=simulate(Gtmp{k},u,'MC',0,'Ts',opt.Ts);
               ymatMC(k,:,:)=ztmp.y;
           end
        else
           ymatMC=[];
        end
    else
        ymatMC=[];
    end
    z=sig(ymat,t,umat,[],ymatMC);
    z.fs=fs;
    xf=xi;    % Final state
else
    %-------------------------------
    %---continuous time model-------
    %-------------------------------
    if ~isnan(u.fs)
        error('LTF.SIMULATE: Signal must be continuous time when model is continuous time')
    end
    umat=u.y;
    [N,nudata]=size(umat);
    if nudata~=nu
        error('LTF.SIMULATE: Number of inputs in SIG object not compatible with model')
    end
    if opt.Ts==0 % Compute default value
        tau=timeconstant(a,fs);  % Dominating pole determines time grid
        Ts=tau/200; % 200 times time constant of system fast enough sampling
        Ts=min([Ts (u.t(end)-u.t(1))/200]);
    else
        Ts=opt.Ts;
    end
    % Common code with ltf.simulate here
    [z,xf]=simulate(lss(G),u,varargin{:});
    z.fs=u.fs;
end

z=inherit(z,G,'Simulation');
if ~isempty(z.name)
    z.name=['Simulation of ',z.name];
end
if nargout==0
    plot(z)
end
end %function


% =================================
% -------Utility functions---------
% =================================

function [b,a,fs]=unpack(s)
%UNPACK fills out b and a polynomials to same orders
%   [b,a,fs]=unpack(s);

[na,nb,nk,nu,ny]=size(s);
fs=s.fs;
b=s.b;
a=s.a;
if nk>0
    for i=1:nu
        for j=1:ny
            bb(j,:,i)=[zeros(1,nk) b(j,:,i)];
        end
    end
    b=bb;
elseif nk<0
    a=[zeros(1,-nk) a];
end
na=length(a);
nb=size(b,2);
if nb>0
    if na>nb
        for i=1:nu
            for j=1:ny
                bbb(j,:,i)=[b(j,:,i) zeros(1,na-nb)];
            end
        end
        b=bbb;
    elseif na<nb
        a=[a zeros(1,nb-na)];
    end
end
end %function


function [na,nb,nk,nu,ny]=size(s,dim)
%SIZE returns the sizes nn=[na,nb,nk,nu,ny]
%   [na,nb,nk,nu,ny]=size(s)
%   Special cases:
%   n=size(s) gives n=[ny nu]. Useful in for instance s+ones(size(s))
%   ny=size(s,1);
%   nu=size(s,2);

nn=s.nn;
na=nn(1); nb=nn(2); nk=nn(3); nu=nn(4); ny=nn(5);
if nargin==2
    if dim==1; nx=ny; end
    if dim==2; nx=nu; end
end
if nargout==1,
    na=[ny nu];
elseif nargout==0,
   disp(['na = ',num2str(na)])
   disp(['nb = ',num2str(nb)])
   disp(['nk = ',num2str(nk)])
   disp(['nu = ',num2str(nu)])
   disp(['ny = ',num2str(ny)])
end
end %function


function texcode=tex(s,varargin)
%TEX creates latex code for LTF objects
%   texcode=tex(s,Property1,Value1,...)
%   The output is the produced texcode that can be pasted into
%   any latex document. Alternatively, the code is put into a file,
%   which is input by reference into the document.
%   Filters (freeware or shareware) for other word processors are available:
%   - texpoint: freeware for Powerpoint
%   - LaImport: FrameMaker
%   - tex2word: Word
%   - latex2rtf: RTF documents
%   - tex4ht: HTML or XML hypertext documents
%
%   Property   Value/{Default}  Description
%   --------------------------------------------------------------
%   format     {'%11.2g'}       Output format, see sprintf for options
%   filename   {''}             Name of the .tex file (none for '')
%   env        {'eqnarray*'}    Tex environment, '' means no env
%
%   Examples:
%     m=rand(ltf([2 1 1]))
%     tex(m)
%   See also: textable, texmatrix, lss.tex

%##


opt=struct('filename','','decimals',1,'env','eqnarray*','format','%11.2g');
opt=optset(opt,varargin);
if ~isempty(opt.env)
    texcode=sprintf(['\\begin{',opt.env,'}']);
else
    texcode='';
end

a=s.a; b=s.b; fs=s.fs;
[na,nb,nk,nu,ny]=size(s);

if isinf(fs)
    var='s';
else
    var='z';
end
for i=1:nu
    for j=1:ny
        astr=mchar(a,var,opt.format,'\cdot ');
        bstr=mchar(b(j,:,i),var,opt.format,'\cdot ');
        nktmp=nk+size(b,2)-length(a);  % Compensate for q and s formalism
        if nktmp>0
            astr=[var,'^',num2str(nktmp),' (',astr,')'];
        elseif nktmp<0
            bstr=[var,'^',num2str(-nktmp),' (',bstr,')'];
        end
        if ny>1
            ystr=[' Y_',int2str(j),'(',var,') &=& '];
        else
            ystr=[' Y(',var,') &=& '];
        end
        if nu>1
            ustr=[' U_',int2str(i),'(',var,') '];
        else
            ustr=[' U(',var,') '];
        end
        if i<nu | j<ny; ustr=[ustr,' \\']; end
        texcode=strvcat(texcode,[ystr,'\frac{',bstr,'}{',astr,'} ',ustr]);
    end
end
if ~isempty(opt.env)
    texcode=strvcat(texcode,sprintf(['\\end{',opt.env,'}']));
end
if ~isempty(opt.filename)
  eval(['fid=fopen(''',opt.filename,'.tex'',''w'')']);
  for j=1:size(texcode,1);
      rowj=texcode(j,:);
      rowj=strrep(rowj,'\','\\');
      fprintf(fid,[rowj,' \n']);
  end
  fclose(fid);
end
end %function

function disp(m)
%DISP gives an ASCII-formatted printout of the LTF object
%   disp(sys)  or  sys
format='%11.2g';

%disp('')
na=m.nn(1); nb=m.nn(2); nk=m.nn(3); nu=m.nn(4); ny=m.nn(5);
fs=m.fs;
b=m.b; a=m.a;
if isempty(b) & isempty(a)
    disp(['Unspecified LTF model with na=',num2str(na),', nb= ',num2str(nb),...
          ', nk= ',num2str(nk),' and with ',...
          num2str(nu),' inputs and ',...
          num2str(ny),' outputs. '])
    return
end

if isnan(m.fs)
    %disp('Continuous time transfer function.')
else
    %disp(['Discrete time transfer function.'])
    fsstr=num2str(m.fs);
    %disp(['Sampling frequency: ',fsstr])
end

if isnan(m.fs)
    var='s';
else
    var='z';
end
ind=find(abs(a)<1e-14); a(ind)=zeros(size(ind));
ind=find(abs(b)<1e-14); b(ind)=zeros(size(ind));
for j=1:ny
    for i=1:nu
        %systmp=minreal(this(j,i));
        astr=mchar(a,var,format);
        bstr=mchar(b(j,:,i),var,format);
        if ny>1
            ystr=[' Y',int2str(j),'(',var,') = '];
        else
            ystr=[' Y(',var,') = '];
        end
        if nu>1
            ustr=[' U',int2str(i),'(',var,')'];
        else
            ustr=[' U(',var,')'];
        end

        ystr=[blanks(length(ystr));ystr;blanks(length(ystr))];
        ustr=[blanks(length(ustr));ustr;blanks(length(ustr))];
        nktmp=nk+size(b,2)-length(a);  % Compensate for q and s formalism
        if nktmp>0
            if nktmp>1
               exponent=['^',num2str(nktmp)];
            else
               exponent='';
            end
            if strcmp(astr,'1');
                astr=[var,exponent];
            else
                astr=[var,exponent,'(',astr,')'];
            end
        elseif nktmp<0
            if nktmp<-1
               exponent=['^',num2str(-nktmp)];
            else
               exponent='';
            end
            if ~strcmp(bstr,'0')
                if strcmp(bstr,'1');
                    bstr=[var,exponent];
                else
                    bstr=[var,exponent,'(',bstr,')'];
                end
            end
        end
        nastr=length(astr);
        nbstr=length(bstr);
        if strcmp(astr,'1');
            if size(b,2)>1
                G=[blanks(nbstr+2);'(',bstr,')';blanks(nbstr+2);];
            else
                G=[blanks(nbstr);bstr;blanks(nbstr);];
            end
        else
            nf=max([nastr,nbstr]);
            if nastr>nbstr
                bstr=[blanks(floor((nf-nbstr)/2)) bstr blanks(ceil((nf-nbstr)/2))];
            else
                astr=[blanks(floor((nf-nastr)/2)) astr blanks(ceil((nf-nastr)/2))];
            end
            frac=char(45*ones(1,nf));
            G=strvcat(bstr,frac,astr);
        end
        disp('')
        disp([ystr G ustr])
    end
end
disp('')
end %function


function info(s)
%INFO displays the name and description.
disp(['Name:    ',s.name])
if ~isempty(s.desc)
    disp(['Description: '])
    disp(s.desc)
end
if ~isempty(s.ulabel)
    disp('Inputs:')
    for i=1:length(s.ulabel)
        disp(['   u',num2str(i),':  ',s.ulabel{i}])
    end
end
if ~isempty(s.ylabel)
    disp('Outputs:')
    for i=1:length(s.ylabel)
        disp(['   y',num2str(i),':  ',s.ylabel{i}])
    end
end
if ~isempty(s.xlabel)
    disp('States:')
    for i=1:length(s.xlabel)
        disp(['   x',num2str(i),':  ',s.xlabel{i}])
    end
end
end %function


% =================================
% -------Help functions---------
% =================================

function [s1,s2]=ltfcheck(s1,s2)
%LTFCHECK
%  Make sure that both arguments are LTF objects
%  LTI objects are converted to LTF
%  PDFCLASS objects are converted to scalar static systems

if isnumeric(s1)       % Static ltf
   % Do nothing
elseif isa(s1,'ltf')
   %OK
elseif isa(s1,'lti')   % Convert
   s1=ltf(s1);
elseif isa(s1,'pdfclass');  % Uncertain static ltf
   bMC=rand(s1,s1.MC,1);  % Take random numbers
   for i=1:s1.MC
      sysMC{i}=ltf(bMC(i),1,s2.fs);
   end
   s1=ltf(E(s1),1,s2.fs);  % Nominal system
   s1.sysMC=sysMC;
else
   error('LTF: Unsupported object type')
end

if isnumeric(s2)       % Static ltf
   % Do nothing
elseif isa(s2,'ltf')
   %OK
elseif isa(s2,'lti')   % Convert
   s2=ltf(s2);
elseif isa(s2,'pdfclass');  % Uncertain static ltf
   bMC=rand(s2,s2.MC,1);  % Take random numbers
   for i=1:s2.MC
      sysMC{i}=ltf(bMC(i),1,s1.fs);
   end
   s2=ltf(E(s2),1,s1.fs);  % Nominal system
   s2.sysMC=sysMC;
end

end %function


function fs=fscheck(s1,s2)
% Same sampling frequency
if isnan(s1.fs) & isnan(s2.fs)   % Continuous time
    fs=NaN;
elseif isnan(s1.fs) | isnan(s2.fs)
    error('LTF: cannot mix continuous and discrete time models')
elseif s1.fs==s2.fs
    fs=s1.fs;
else
    error('LTF: Multi-rate LTI operations are not possible')
end


end %function


function sysMC=mceval1(op,s,varargin)
%MCEVAL1 utility function to compute MC data
if isa(s,'ltf') & ~isempty(s.sysMC)
    for i=1:length(s.sysMC)
        sysMC{i}=feval(op,s.sysMC{i},varargin{:});
    end
else
    sysMC=[];
end

end %function


function sysMC=mceval2(op,s1,s2)
%MCEVAL2 utility function to compute MC data
if isa(s1,'lti') & ~isempty(s1.sysMC)
    MC1=length(s1.sysMC);
else
    MC1=0;

end
if isa(s2,'lti') & ~isempty(s2.sysMC)
    MC2=length(s2.sysMC);
else
    MC2=0;
end
MC=max([MC1 MC2]);
if MC>0
    if isa(s1,'lti')
       s1tmp=rand(s1,MC);
    else
       s1tmp{1:MC}=s1;
    end
    if isa(s2,'lti')
       s2tmp=rand(s2,MC);
    else
       s2tmp{1:MC}=s2;
    end
    for i=1:MC
        sysMC{i}=feval(op,s1tmp{i},s2tmp{i});
    end
else
    sysMC=[];
end
end %function

end %methods
end %ltf
