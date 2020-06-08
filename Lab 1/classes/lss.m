classdef lss < ltimod
%LSS State-Space implementation of Linear Time Invariant (LTI) systems
%   Causal MIMO deterministic and stochastic systems supported.
%
%   lss.A, lss.B         State dynamics Ax+Bu
%   lss.C, lss.D         Measurement relation Cx+Du
%   lss.Q, lss.R, lss.S  Covariances Cov(v),Cov(e), Cov(v,e)
%   lss.fs               Sampling frequency
%
%   help lss.lss      gives help for the constructor
%   methods lss       lists all methods for this class

%   Copyright Fredrik Gustafsson, Sigmoid AB
%   $ Revision: 28-Oct-2019 $

properties (SetAccess = public)
  A; B; C; D; Q; R; S; nn;
  fs; name; xlabel; ulabel; ylabel; tlabel; marker; markerlabel; desc;
  MC; pe; pv;
  sysMC;  % Monte Carlo samples
end

methods
% =================================
% -----------Constructor-----------
% =================================

function m=lss(varargin)
%LSS Construct a State-Space object
%   Model:
%      x(t+1) = A x(t) + B u(t) + v(t)
%       y(t)  = C x(t) + D u(t) + e(t)
%       Q=E(vv'), R=E(ee') and S=E(ve')
%
%   The orders are denoted nn=[nx nu nv ny] in 'order of appearance'
%   which happens to be the alphabetical order.
%   The noise distributions can be set with the fields pv and pe, resp.
%   For instance, s.pe=udist(-1,1); gives uniform noise on output
%
%   Assignment for deterministic input-output dynamics:
%     LSS(A,B,C,D)      Continuous systems (fs=NaN)
%     LSS(A,B,C,D,fs)   Discrete systems
%
%   Assignment for stochastic systems:
%     LSS(A,B,C,D,Q,R,S,fs) or
%     LSS(A,B,C,D,Q,R,fs)   if S=0
%
%   For noise models, use
%     LSS(A,[],C,[],Q,R,S,fs)
%   For static noise-free models, use
%     LSS([],[],[],D,fs)
%   Empty structure:
%     LSS([nx,nu,nv,ny])
%     LSS([nx,nu,nv,ny],fs)
%
%   Special constructions:
%     LSS([nx,nu,nv,ny]) defines state space structure for rand and estimate
%     LSS('unit')  defines the unit system y=u
%     LSS('delay') defines the unit delay system y(t)=u(t-1)
%     LSS('sum')   defines the summator (integrator approx.) y(t)=y(t-1)+u(t)
%     LSS('int')   defines the integrator G(s)=1/s
%
%   Examples:
%     mc=lss([1 1;0 1],[0;1],[1 0],0)
%     md=lss([1 1;0 1],[0;1],[1 0],0,2)
%     ms=lss([1 1;0 1],[0;1],[1 0],0,[0 0;0 1],1,2)
%   See also: LTF

%#
%#


A=[]; B=[]; C=[]; D=[]; Q=[]; R=[]; S=[];
fs=NaN;
name='';  xlabel=[]; ulabel=[]; ylabel=[]; desc=[];
nn=[0 0 0 0];
MC=0;
pv=[]; pe=[];
sysMC=[];

if nargin<1
%    error('LSS contructor expects an input argument: see help lss')
elseif isa(varargin{1},'lss')|isa(varargin{1},'nl')|isa(varargin{1},'ltf')||isa(varargin{1},'arx')
    s=varargin{1};
    if nargin>1
       MC=varargin{2};
    else
       MC=30;
    end
    if isa(s,'lss')  % Do nothing
        s2=s;
    elseif isa(s,'nl')  % NL object
        if nargin<2
           error('Conversion from NL to LSS requires a linearization point')
        end
        s2=nl2lss(s,varargin{2:end});
    elseif isa(s,'ltf')  % Convert from LTF object
        s2=ltf2lss(s,varargin{2:end});
        % Compute MC samples
        if ~isempty(s.sysMC)
            MC=length(s.sysMC);
            srand=rand(s,MC);            %MC samples of s...
            if length(srand)==MC % Empty cell array may be returned
               for k=1:MC
                  sysMC{k}=lss(srand{k});    %are converted to lss
               end
            else
               MC=0;
            end
        else
            MC=0;
        end
    elseif isa(s,'arx')  % Convert from ARX object
        s2=arx2lss(s);
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
                sysMC{k}=lss(srand{k});     %are converted to lss
            end
        else
            MC=0;
        end
    end
    A=s2.A; B=s2.B; C=s2.C; D=s2.D; Q=s2.Q; R=s2.R; S=s2.S;
    fs=s2.fs; nn=s2.nn;
    name=s2.name; ulabel=s2.ulabel; ylabel=s2.ylabel; desc=s2.desc;
elseif nargin==1 | nargin==2
    if ~isstr(varargin{1})
        if nargin==1;
            fs=NaN;
        else
            fs=varargin{2};
            if ~isnumericscalar(fs) & ~isnan(fs)
                error('LSS constructor: fs expected as second input argument')
            end
            if fs==0 | isinf(fs)
               fs=NaN;
            elseif fs<0
               error('LSS constructor: fs>0')
            end
        end
        nn=varargin{1};
        if ~isnumeric(nn) | length(nn)~=4
            error('LSS constructor: Structure index must be of length(4), nn=[nx nu nv ny]')
        else
            nx=nn(1); nu=nn(2); nv=nn(3);  ny=nn(4);
            if nv>0; ne=ny; end;
        end
        if ny==0
            error('LSS constructor: ny>0 required')
        end
        if any(nn<0)
            error('LSS constructor: orders nn must be positive')
        end
        if any(round(nn)~=nn)
            error('LSS constructor: orders nn must be integers')
        end
%        error('LSS constructor: String input required when single argument used to LSS')
    elseif isstr(varargin{1})
        if strcmpi(varargin{1},'unit')
            s=lss([],[],[],1);
        elseif strcmpi(varargin{1},'sum')
            s=lss(1,1,1,1,1);
        elseif strcmpi(varargin{1},'delay')
            s=lss(0,1,1,0,1);
        elseif strcmpi(varargin{1},'int')
            s=lss(0,1,1,0,NaN);
        else
            error('Unknown string option to LSS')
        end
        if nargin>1
            error(sprintf('''%s'' may not be followed by any arguments.', ...
                  varargin{1}))
        end
        A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs;
        [ny,nu] = size(D);
        nx=size(A,1);
        nv=rank(Q);
    end
    nn=[nx nu nv ny];
else  % ss(A,B,C,D,...)
    A=varargin{1};
    B=varargin{2};
    C=varargin{3};
    if nargin>3
       D=varargin{4};
    else
       D=[];
    end
    if nargin==4
        fs = NaN;
    end
    [Q,R,S,fs]=lssdatacheck(lss,A,B,C,D,varargin{5:end});
    ny=size(C,1);
    nu=size(B,2);
    if ~isempty(D)
       [ny,nu] = size(D);
    else
       D=zeros(ny,nu);
    end
    if isempty(A)
       nx=0;
    else
       nx=size(A,1);
    end
    nv=rank(Q);
    nn=[nx nu nv ny];
end
if size(B,2)==0, B=zeros(nn(1),0); end
if size(D,2)==0, D=zeros(nn(4),0); end

m.A=A; m.B=B; m.C=C; m.D=D; m.Q=Q; m.R=R; m.S=S;
m.fs=fs; m.nn=nn;
m.name=name; m.ulabel=ulabel; m.ylabel=ylabel; m.desc=desc;
end


function out=fieldread(arg)
out=eval(arg);
end

function fieldwrite(arg1,arg2)
switch arg1
    case {'A', 'B', 'C', 'D', 'Q', 'R', 'S', 'nn'}
        if ~isnumeric(arg2)
            error('Field value must be numeric.')
        end
end
if strcmp(arg1,'sysMC')
      if isempty(arg2)
         sysMC=[];
         MC=0;
      elseif ~isa(arg2,'cell')
         error('sysMC must be a cell array')
      elseif ~isa(arg2{1},'lss')
         error('sysMC must be a cell array of LSS objects')
      elseif isequal(arg2{1}.nn,nn)  % Check sysMC argument
         sysMC=arg2;
         MC=length(sysMC);
      else
         error('Argument of sysMC does not match the LSS dimensions')
      end
elseif ~all(size(arg2)==size(eval(arg1)))
   error('Cannot change size of protected fields in LSS objects')
else
   eval([arg1,'=',mat2str(arg2),';'])
end
end

function sys=arrayread(j,i)
%ARRAYREAD used to pick out sub-systems by indexing
%   sys=arrayread(j,i)
%   j is the row index/indices, corresponding to the outputs
%   i is the column index/indices, corresponding to the inputs
%   Example: s(2,3) gives the SISO system from input 3 to output 2

% see csl file
end

% =================================
% ----------Concatenations---------
% =================================

function sys=ctranspose(s)
%CTRANSPOSE reverse the inputs with outputs for MIMO systems
%   sys=ctranspose(s)  or sys=s'
%   ctranspose is the same as transpose
sys=transpose(s);
end

function sys=transpose(s)
%TRANSPOSE reverse the inputs with outputs for MIMO systems
%   sys=transpose(s)  or sys=s.'
%
%   For state space models, this corresponds to
%   (A,B,C,D) -> (A',C',B',D')

sys=lss(s.A',s.C',s.B',s.D',s.fs);
sys.sysMC=mceval1('transpose',s);
end

function sys=horzcat(varargin)
%HORZCAT performs horizontal concatenation
%   sys=horzcat(s1,s2,...)  or  [s1 s2 ...]
%   This can be used to create a MISO system from SISO systems
%   Reguirement: All systems must have the same number of outputs ny1=ny2=...

if nargin<1;
    error('No system input argument')
elseif nargin<2;
    sys=varargin{1};
else % Concatenate the first two systems
    s1=varargin{1};
    s2=varargin{2};
    if isa(s1,'lss') | isa(s2,'lss')
        [s1,s2]=lsscheck(s1,s2);
        [nx1,nu1,nv1,ny1]=size(s1);
        [nx2,nu2,nv2,ny2]=size(s2);
        if ny1~=ny2
            error(['Incorrect dimensions in horzcat. ny1=ny2 required, but ny1=',...
                    num2str(ny1),' and ny2=',num2str(ny2)])
        end
        A=[s1.A zeros(nx1,nx2); zeros(nx2,nx1) s2.A];
        B=[s1.B zeros(nx1,nu2); zeros(nx2,nu1) s2.B];
        C=[s1.C s2.C];
        D=[s1.D s2.D];
        Q=[s1.Q zeros(nx1,nx2); zeros(nx2,nx1) s2.Q];
        R=s1.R+s2.R;
        S=[s1.S;s2.S];
        sys=lss(A,B,C,D,Q,R,S,s1.fs);
    end
end
if nargin>2 % Recursion for the remaining elements
   sys=horzcat(sys,varargin{3:end});
end
sys.name=['Concatenation of ',s1.name,' and ',s2.name];
sys.desc=char(s1.desc,s2.desc);
if ~isempty(s1.ulabel) | ~isempty(s2.ulabel);
   sys.ulabel={s1.ulabel{1:nu1},s2.ulabel{1:nu2}};
end
if ~isempty(s1.ylabel)
   sys.ylabel=s1.ylabel;  % Same outputs for vertcat
end
if ~isempty(s1.xlabel) | ~isempty(s2.xlabel);
   sys.xlabel={s1.xlabel{1:nx1},s2.xlabel{1:nx2}};
end
sys.sysMC=mceval1('horzcat',varargin{:});
end


function sys=vertcat(varargin)
%VERTCAT performs vertical concatenation
%   sys=vertcat(s1,s2,...) or [s1;s2;...]
%   This can be used to create a SIMO system from SISO systems
%   Reguirement: All systems must have the same number of inputs nu1=nu2=...

if nargin<1;
    error('No system input argument')
elseif nargin<2;
    sys=varargin{1};
else % Concatenate the first two systems
    s1=varargin{1};
    s2=varargin{2};
    if isa(s1,'lss') | isa(s2,'lss')
        [s1,s2]=lsscheck(s1,s2);
        [nx1,nu1,nv1,ny1]=size(s1);
        [nx2,nu2,nv2,ny2]=size(s2);
        if nu1~=nu2
            error(['Incorrect dimensions in horzcat. nu1=nu2 required, but nu1=',...
                    num2str(nu1),' and nu2=',num2str(nu2)])
        end
        A=[s1.A zeros(nx1,nx2); zeros(nx2,nx1) s2.A];
        C=[s1.C zeros(ny1,nx2); zeros(ny2,nx1) s2.C];  % different outputs
        B=[s1.B; s2.B];  % same input
        D=[s1.D; s2.D];
        Q=[s1.Q zeros(nx1,nx2); zeros(nx2,nx1) s2.Q]; % different process noises
        R=[s1.R zeros(ny1,ny2); zeros(ny2,ny1) s2.R]; % different measurement noises
        S=[s1.S zeros(nx1,ny2); zeros(nx2,ny1) s2.S];
        sys=lss(A,B,C,D,Q,R,S,s1.fs);
    end
end
if nargin>2 % Recursion for the remaining elements
   sys=vertcat(sys,varargin{3:end});
end
sys.name=['Concatenation of ',s1.name,' and ',s2.name];
sys.desc=char(s1.desc,s2.desc);
if ~isempty(s1.ylabel) | ~isempty(s2.ylabel);
   sys.ylabel={s1.ylabel{1:ny1},s2.ylabel{1:ny2}};
end
if ~isempty(s1.ulabel)
   sys.ulabel=s1.ulabel;  % Same inputs for vertcat
end
if ~isempty(s1.xlabel) | ~isempty(s2.xlabel);
   sys.xlabel={s1.xlabel{1:nx1},s2.xlabel{1:nx2}};
end
sys.sysMC=mceval1('vertcat',varargin{:});
end


% =================================
% -----Over-loaded operators-------
% =================================



%------- Functions of one system-------------

function sys=uplus(sys)
%UPLUS unitary plus
%   sys=uplus(sys) or sys=+sys

% do nothing
end

function sys=uminus(sys)
%UMINUS  unitary minus
%   sys=uminus(sys) or sys=-sys

sys.B = -sys.B;
sys.D = -sys.D;
sys.sysMC=mceval1('uminus',sys);
end

function sys=minus(s1,s2)
%MINUS parallel connection with difference at output
%   sys=minus(s1,s2) or sys=s1-s2
%   Requirements: nu1=nu2, ny1=ny2

sys = s1+(-s2);
sys.sysMC=mceval2('minus',s1,s2);
end

function sys=power(s,n)
%POWER is repeated multiplication
%   sys=power(s,n)  or  s.^n
%   calls s^n
if nargin<2; error('LSS.POWER needs two input arguments'), end
sys=s^n;
sys.sysMC=mceval1('power',sys,n);
end

function sys=mpower(s,n)
%MPOWER is repeated multiplication
%   sys=mpower(s,n)  or  s^n

if nargin<2; error('LSS.MPOWER needs two input arguments'), end
if round(n)~=n
    error('Power must be an integer')
end
if ~isreal(n)
    error('Power must be a real number')
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
sys.sysMC=mceval1('mpower',sys,n);
end

%------- Functions of two systems-------------

function ans=eq(s1,s2)
%EQ compares two systems for equality
%   ans=eq(s1,s2)  or  ans=(s1==s2)

if isa(s1,'lti') & ~isa(s1,'lss')
   s1=lss(s1);
end
if isa(s2,'lti') & ~isa(s2,'lss')
   s2=lss(s2);
end

if isa(s1,'lss') & isa(s2,'lss')
    % Compare nx impulse response coefficients
    comp=(norm(s1.D-s2.D)<1e-12);
    nx=s1.nn(1);
    for i=1:nx;
        comp= comp & (norm(s1.C*s1.A^i*s1.B - s2.C*s2.A^i*s2.B)<1e-12);
    end
    % compare covariance matrices
    if s1.nn(3)>0 | s2.nn(3)>0
       q=(norm(s1.Q - s2.Q)<1e-12);
       r=(norm(s1.R - s2.R)<1e-12);
       s=(norm(s1.S - s2.S)<1e-12);
    else
       q=1; r=1; s=1;
    end
    ans=false;
    if comp & q & s & r
        ans=true;
    end
else
    error('Both s1 and s2 must be LSS objects for EQ')
end
end

function b=isequal(s1,s2,conf);
%ISEQUAL performs a hypothesis test with confidence level conf
%   b=isequal(s1,s2,conf);
%   Note: Only the input-output part is tested (A,B,C,D)
%   b is true or false
%   conf is the confidence level (default 0.99)
%
%   See also: eq (non-stochastic version)

if nargin<3, conf=0.99; end

if any(s1.nn~=s2.nn)
    b=false;
else
    m1=arx(ltf(s1));
    m2=arx(ltf(s2));
    b=isequal(m1,m2,conf);
end
end

function sys=mrdivide(s1,s2)
%MRDIVIDE computes the right inverse of a system such that s*sys=I
%   sys=s1*inv(s2) is used
sys=s1*inv(s2);
sys.sysMC=mceval2('mrdivide',s1,s2);
end

function sys=mldivide(s1,s2)
%MLDIVIDE computes the left inverse of a system such that s*sys=I
%   sys=inv(s1)*s2 is used
sys=inv(s1)*s2;
sys.sysMC=mceval2('mldivide',s1,s2);
end

function sys=plus(s1,s2)
%PLUS Add two models together in a parallel connection
%   sys=plus(s1,s2)  or sys=s1+s2
%   Requirements: nu1=nu2, ny1=ny2

if isnumericscalar(s1)
    s2.D = s2.D+s1;
    sys = s2;
    return
elseif isnumericscalar(s2)
    s1.D = s1.D+s2;
    sys = s1;
    return
else
    [s1,s2]=lsscheck(s1,s2);
    [nx1,nu1,nv1,ny1]=size(s1);
    [nx2,nu2,nv2,ny2]=size(s2);
    if  ny1~=ny2
        error('Systems s1 and s2 are not compatible for parallel connection (+): ny1~=ny2')
    end
    if  nu1~=nu2
        error('Systems s1 and s2 are not compatible for parallel connection (+): nu1~=nu2')
    end
    a = [s1.A zeros(nx1,nx2); zeros(nx2,nx1) s2.A];
    sys =lss(a, [s1.B; s2.B], [s1.C s2.C], s1.D+s2.D, s1.fs);
end
sys.sysMC=mceval2('plus',s1,s2);
end

function sys=mtimes(s1,s2)
%MTIMES computes series connection of systems
%   sys=mtimes(s1,s2)  or   sys = s1*s2
%   IMplements the relation y(t)=s1*s2*u(t)
%   Requirements: ny2=nu1

if isnumeric(s1)
    % y=c*G*u scales output
    A=s2.A;
    B=s2.B;
    C = s1*s2.C;
    D = s1*s2.D;
    fs=s2.fs;
elseif isnumeric(s2)
    % y=G*c*u scales input
    A=s1.A;
    B = s1.B*s2;
    C=s1.C;
    D = s1.D*s2;
    fs=s1.fs;
else
    [s1,s2]=lsscheck(s1,s2);
    [nx1,nu1,nv1,ny1]=size(s1);
    [nx2,nu2,nv2,ny2]=size(s2);
    fs=s1.fs;
    if  ny2~=nu1
        error(['In series connection (*), ny2~=nu1, but ny2=',num2str(ny2),' nu1=',num2str(nu1)])
    end

%    a = [s1.A zeros(nx1,nx2); s2.B*s1.C s2.A];
%    b = [s1.B; s2.B*s1.D];
%    c = [s2.D*s1.C s2.C];
%    d = s2.D*s1.D;
    A = [s1.A s1.B*s2.C; zeros(nx2,nx1) s2.A];
    B = [s1.B*s2.D; s2.B];
    C = [s1.C s1.D*s2.C];
    D = s1.D*s2.D;
end
sys=lss(A,B,C,D,fs);
sys.sysMC=mceval2('mtimes',s1,s2);
end

function sys=inv(s)
%INV computes the right inverse of a system such that s*sys=I
%   sys=inv(s)
%   Requirements: ny<=nu and rank(D)=ny
[nx,nu,nv,ny]=size(s);
if nv>0
    error('No inverse for stochastic systems')
end
A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs;
Dinv=D\eye(ny);
A1=A-B*Dinv*C;
B1=B*Dinv;
C1=-Dinv*C;
D1=Dinv;
sys=lss(A1,B1,C1,D1,fs);
end

% =================================
% ----------LTI functions----------
% =================================

function sys=transform(s,T)
%TRANSFORM computes the transformed LSS model with state T*x
%   sys=transform(s)
%
%   Requirements: size(T)=[nx nx]

[nx,nu,nv,ny]=size(s);
if ~isequal(size(T),[nx nx])
    error('LSS.transform: T must be square with size nx')
end
if rank(T)<nx
    error('LSS.transform: T must have full rank')
end
Tinv=inv(T);
if nv==0
   sys=lss(T*s.A*Tinv,T*s.B,s.C*Tinv,s.D,s.fs);
else
   sys=lss(T*s.A*Tinv,T*s.B,s.C*Tinv,s.D,T*s.Q*T',s.R,T*s.S,s.fs);
end
end

function sys=feedback(s1,s2)
%FEEDBACK computes the feedback connection of two systems
%   sys=feedback(s1,s2)
%   System s1 in forward loop and s2 in feedback loop.
%   Requirements: nu1=ny2, nu2=ny1

[s1,s2]=lsscheck(s1,s2);
[nx1,nu1,nv1,ny1]=size(s1);
[nx2,nu2,nv2,ny2]=size(s2);

if  ny1~=nu2
    error('Systems s1 and s2 are not compatible for feedback connection: ny1~=nu2')
end
if  ny2~=nu1
    error('Systems s1 and s2 are not compatible for feedback connection: ny2~=nu1')
end
if nx2==0 % exception for feedback path
    x = eye(ny2)+s2.D*s1.D;
    a = [s1.A-s1.B*(x\s2.D)*s1.C];
    b = [s1.B/x];
    c = [(eye(ny1)-s1.D*(x\s2.D))*s1.C ];
    d = s1.D/x;
elseif all(all(s1.D==0)) & all(all(s2.D==0))
    a = [s1.A -s1.B*s2.C; s2.B*s1.C s1.A];
    b = [s1.B; zeros(nx2,nu2)];
    c = [s1.C zeros(ny2,nx2)];
    d = zeros(ny1,ny1);
else
    x = eye(ny2)+s2.D*s1.D;
    a = [s1.A-s1.B*(x\s2.D)*s1.C,             -s1.B*(x\s2.C);
         s2.B*(eye(ny1)-s1.D*(x\s2.D))*s1.C,  s2.A-s2.B*s1.D*(x\s2.C)];
    b = [s1.B/x; s2.B*s1.D/x];
    c = [(eye(ny1)-s1.D*(x\s2.D))*s1.C -s1.D*(x\s2.C)];
    d = s1.D/x;
end
sys =lss(a,b,c,d,s1.fs);
sys.sysMC=mceval2('feedback',s1,s2);
end

function sys=diag(varargin)
%DIAG appends independent models into one using append recursively
%   sys=diag(s1,s2,s3,...)  % Observe, no brackets
%   sys=diag(s)  with one input argument extract the diagonal of the
%   transfer function matrix s

if nargin==1;
   s=varargin{1};
   n=min(s.nn([2 4]));
   sys=s(1,1);
   for k=2:n
      sys=diag(sys,s(k,k));
   end
   sys.sysMC=mceval1('diag',varargin{1});
elseif nargin>1;
   sys=append(varargin{1},varargin{2});
   sys.sysMC=mceval2('append',varargin{1},varargin{2});
end
if nargin>2  % Recursion
   sys=diag(sys,varargin{3:end});
   sys.sysMC=mceval2('diag',sys,varargin{3:end});
end
end

function sys=append(s1,s2)
%APPEND packs two independent models into one.
%   sys = append(sys1,sys2)
%   Requirements: none

[s1,s2]=lsscheck(s1,s2);
[nx1,nu1,nv1,ny1]=size(s1);
[nx2,nu2,nv2,ny2]=size(s2);

a = [s1.A zeros(nx1,nx2); zeros(nx2,nx1) s2.A];
b = [s1.B zeros(nx1,nu2); zeros(nx2,nu1) s2.B];
c = [s1.C zeros(ny2,nx1); zeros(ny1,nx2) s2.C];
d = [s1.D zeros(ny1,nu2); zeros(ny2,nu1) s2.D];
sys =lss(a,b,c,d,s1.fs);
sys.sysMC=mceval2('append',s1,s2);
end

function sys=downsample(s,n)
%DOWNSAMPLE downsamples a discrete time LSS object a factor n
%   sys = downsample(s)
%   The new system has nu*n inputs and ny*n outputs
%   Use the corresponding function for the SIG object
%   Requirements: none

A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S;
ylabel={s.ylabel};
ulabel={s.ulabel};

[U,V]=svd(Q);
nv=length(find(abs(diag(V))>1e-12*V(1,1)));
G=U(:,1:nv);
Q=V(1:nv,1:nv);

Ad=A;
Bd=B;
Qtmp=G*Q*G';
Qd=Qtmp;
Cd=C;
Sd=zeros(size(C));
ydlabel=ylabel;
udlabel=ulabel;
for k=2:n
   Ak=A^(k-1);
   Qtmp=Ak*Qtmp*Ak';
   Qd=Qd+Qtmp;
   Cd=[Cd; C*Ak];
   Bd=[Bd Ak*B];
   Sd=[Sd; C*A^(k-2)*Q*A'];
ydlabel,ylabel
   ydlabel={ydlabel{:},ylabel{:}};
   udlabel={udlabel{:},ulabel{:}};
end
Rd=kron(eye(n),R);
Dd=kron(eye(n),D);
size(Sd)
sys =lss(Ad,Bd,Cd,Dd,Qd,Rd,Sd',s.fs/n);
sys.xlabel=s.xlabel;
sys.ylabel=ydlabel;
sys.ulabel=udlabel;
sys.desc=s.desc;
sys.name=[s.name,' decimated a factor ',num2str(n)];
sys.sysMC=mceval1('decimate',s);
end

function sys=c2d(s,fs,method)
%C2D converts continuous time system to discrete time
%   sys=c2d(s,fs,method)
%
%   fs {1} is the sampling frequency
%   method is a string describing the assumption on intersample behaviour
%   {'ZOH'}    zero order hold, piece-wise constant input assumed
%   'FOH'      first order hold, piecewise linear input assumed
%   'bilinear' s=2/T (z-1)/(z+1)
%
%   References:
%     Astrom and Wittenmark, Computer controlled systems, Prentice-Hall
%     ZOH page 37 (3.6) and (3.7)
%     FOH page 41 (3.18) and (3.19)
%   Examples:
%     mc=lss(ltf(1,[1 1 1]))
%     md=c2d(mc,0.1)

if nargin<3; method='zoh'; end
if nargin<2; fs=1; end

if s.fs>0
    error('LSS object already in discrete time')
end

T=1/fs;
[nx,nu,nv,ny]=size(s);
Ac=s.A; Bc=s.B; Cc=s.C; Dc=s.D;
method=lower(method);
if strcmp(method(1:3),'zoh');
    % ZOH sampling relation: [Ad Bd;zeros(nu,nx) eye(nu)]=expm([ [Ac Bc]*T; zeros(nu,nx+nu)]);
    E=real(expm([ [Ac Bc]*T; zeros(nu,nx+nu)]));
    A=E(1:nx,1:nx);
    C=Cc;
    if nu>0
        B=E(1:nx,nx+1:nx+nu);
        D=Dc;
    else
        B=[]; D=[];
    end
elseif strcmp(method(1:3),'foh');
    E=real(expm([ [Ac Bc]*T zeros(nx,nu); zeros(nu,nx+nu) eye(nu);zeros(nu,nx+2*nu)]));
    E=E(1:nx,:);
    E1=E(:,1:nx); E2=E(:,nx+1:nx+nu); E3=E(:,nx+nu+1:nx+2*nu);
    A=E1;
    C=Cc;
    if nu>0
        B=E2+A*E3-E3;
        D=Dc+Cc*E3;
    else
        B=[]; D=[];
    end
elseif strcmp(method(1:3),'bil');
    if 1  % Use the LTF implementation
      m=lss(c2d(ltf(s),fs,'bil'));
      [A,B,C,D]=getss(m);
    else
       H=inv(eye(size(Ac))-T/2*Ac);
       A=H*(eye(size(Ac))+T/2*Ac);
       C=Cc*H;
       if nu>0
           B=H*T/2*Bc;
           D=Dc+C*B;
       else
           B=[]; D=[];
       end
    end
end
sys =lss(A,B,C,D,s.Q/fs,s.R,s.S/sqrt(fs),fs);
sys=inherit(sys,s,'Discretized model');
sys.sysMC=mceval1('c2d',sys,fs,method);
end

function sys=d2c(s,method)
%D2C converts discrete time system to continuous time
%   sys=d2c(s,method)
%
%   method is a string describing the assumption on intersample behaviour
%   {'ZOH'}    zero order hold, piece-wise constant input assumed
%   'FOH'      first order hold, piece-wise constant input assumed
%   'bilinear' the ltf method d2c is called
%   'ssbilinear' s=2/T (z-1)/(z+1) used in the state space model
%              Both 'bil' and 'ssbil' gives same result but with
%              different state coordinates
%   References:
%     Astrom and Wittenmark, Computer controlled systems, Prentice-Hall
%     ZOH page 39 (3.14) and before
%
%   Example:
%     Gcltf=getfilter(2,0.5,'fs',0);
%     Gctmp=lss(Gcltf)
%     Gd=c2d(Gctmp,1,'zoh')  % Discrete time prototype
%     Gc=d2c(Gd,'zoh');
%     subplot(2,1,1), plot(Gc,Gd,'plottype','plot')
%     subplot(2,1,2), plot(step(Gc,5),step(Gd,5))


if nargin<2; method='zoh'; end
if isnan(s.fs)
    error('LSS object already in continuous time')
end

[nx,nu,nv,ny]=size(s);
Ad=s.A; Bd=s.B; Cd=s.C; Dd=s.D; fs=s.fs;
if isnan(fs) | ~isnumericscalar(fs)
    error('fs must be specified as a scalar')
end
T=1/fs;
method=lower(method);
if strcmp(method(1:3),'zoh');
    % Inverse of ZOH sampling
    % ZOH sampling relation: [Ad Bd;zeros(nu,nx) eye(nu)]=expm([ [Ac Bc]*T; zeros(nu,nx+nu)]);
    E=real(logm([Ad Bd;zeros(nu,nx) eye(nu)]));
    A=fs*E(1:nx,1:nx);
    C=Cd;
    B=fs*E(1:nx,nx+1:nx+nu);
    D=Dd;
elseif strcmp(method(1:3),'foh');
    error('d2c with FOH is not available')
elseif strcmp(method(1:3),'bil');
        % Use the LTF implementation
       m=lss(d2c(ltf(s),'bil'));
       [A,B,C,D]=getss(m);
elseif strcmp(method(1:3),'ssb');
       % use the direct state space transformation
       A=2/T*inv(eye(size(Ad))+Ad)*(Ad-eye(size(Ad)));
       Hinv=(eye(size(A))-T/2*A);
       B=2/T*Hinv*Bd;
       C=Cd*Hinv;
       D=Dd-Cd*Bd;
else
    error('Unknown method')
end
sys =lss(A,B,C,D,s.Q*fs,s.R,sqrt(fs)*s.S,NaN);
sys=inherit(sys,s,'Continuous model computed from discrete model');
sys.sysMC=mceval1('c2d',sys,method);
end

function [sys,T,uind,yind]=minreal(s)
%MINREAL computes a minimal realization
%   [sys,T,uind,yind]=minreal(s)
%
%   A model reduction without information loss is performed
%   using modred.
%   Unused inputs and outputs are removed in a second step.
%   T is the transformation matrix such that xnew=T*x
%   uind and yind are the removed inputs and outputs, respectively.
%
%   Algorithm:
%     1. lss.modred is called with automatic choice of model order,
%     so the trunction is done where all eigenvalues of the balanced
%     A matrix are exactly zero.
%     2. The columns of the B and D matrix are screened for all zero vectors,
%     corresponding to unused inputs, which in such case are removed.
%     The C matrix is screened for zero columns, which are removed.
%
%   Examples:
%     G=lss(diag([-1 -2 -3]),[1;1;0],[1 0 1],0,0)
%     ltf(G)
%     minreal(G)
%     G=lss(diag([-1 -2 -3]),diag([1;1;0]),diag([1 0 1]),zeros(3),0)
%     minreal(G)


% Step 1: loss-free truncatation of balanced realization
[s,T]=modred(s);

% Step 2: Check input and output dependences
A=s.A; B=s.B; C=s.C; D=s.D; Q1=s.Q; R1=s.R; S1=s.S; fs=s.fs; nn=s.nn;
nu=nn(2); ny=nn(4);
uind=[];
for i=1:nu;
    if all(abs(B(:,i))<1e-14) & all(abs(D(:,i))<1e-14)
        uind=[uind i];
    end
end
B(:,uind)=[];
D(:,uind)=[];
yind=[];
for i=1:ny;
    if all(abs(C(i,:))<1e-14) & all(abs(D(i,:))<1e-14)
        yind=[yind i];
    end
end
C(yind,:)=[];
D(yind,:)=[];

% Step 3: Save model
if nn(3)>0 %Stochastic system
    sys=lss(A,B,C,D,T*Q1*inv(T),R1,T*S1,fs);
else
    sys=lss(A,B,C,D,fs);
end
sys=inherit(sys,s,'Minimal realization');
sys.sysMC=mceval1('minreal',sys);
end

function [sys,T]=modred(s,n)
%MODRED computes a reduced order model using a balanced realization
%   [mred,T]=modred(m,n);
%
%   The function first calls balreal to get a balanced realization,
%   where the eigenvalues of the A matrix are sorted in order.
%   This model is then truncated at order n.
%
%   T is the transformation matrix such that xnew=T*x
%
%   n is the order of the filter after model reduction
%   if n is empty or zero, the model order is chosen automatically
%   in which case modred and minreal is the same thing.
%
%
%   Examples:
%     m=rand(lss([8,1,0,1]));
%     mred2=modred(m,2);
%     mredauto=modred(m,0);
%     plot(m,mred2,mredauto);
%
%   See also: gram, balreal

%$ Revision: 28-Oct-2019 $

% Step 1: compute balanced realization
[s,T,S]=balreal(s);
A=s.A; B=s.B; C=s.C; D=s.D; Q1=s.Q; R1=s.R; S1=s.S; fs=s.fs; nn=s.nn;

% Step 2: Analyze eigenvalues in S
if nargin<2 | isempty(n) | n==0  % auto-truncate
    s1=diag(S);
    ind=find(abs(s1/s1(1))<1e-2);
    if isempty(ind)
        n=length(s1);
    else
        n=ind(1)-1;
    end
end

% Step 3: Truncate
Ared=real(A(1:n,1:n));
Bred=real(B(1:n,:));
Cred=real(C(:,1:n));
Dred=real(D);
T=real(T(1:n,:));  % xnew=T*x

% Step 4: Save model
if nn(3)>0 %Stochastic system
    sys=lss(Ared,Bred,Cred,Dred,T*Q1*T',R1,T*S1,fs);
else
    sys=lss(Ared,Bred,Cred,Dred,fs);
end
sys=inherit(sys,s,'Model reduction');
sys.sysMC=mceval1('modred',sys,n);
end

function [sys,T,S]=balreal(s)
%BALREAL computes the balanced realization
%   [sys,T,S]=balreal(s)
%
%   A balanced realization is defined as a state space realization
%   having the same diagonal input and output gramian
%         gram(A,B)=gram(A',C') = diagonal
%   The transfer function is unchanged by the change of state variables
%   xnew=T*x, where the transformation matrix is provided as the second output.
%   S is the controllability and observability matrix of the balanced realization,
%   that can be interpreted as a singular value matrix that is used for
%   model reduction and approximation purposes.
%   Balanced realizations are used for avoiding numerical problems and
%   for model approximation purposes as done in modred.
%   It can also be used to initialize ARMA models from
%   high order AR models in sig2lti.
%
%   A,B,C denote the initial state space model of any realization of a LTF
%   Ab,Bb,Cb denote the state space model of a balanced realization
%   T is the state transformation matrix
%
%   Examples:
%     m=rand(lss([2 1 0 1]))
%     mbal=balreal(m)
%     gram(mbal,'o')
%     gram(mbal,'c')
%
%   See also: gram, modred

A=s.A; B=s.B; C=s.C; D=s.D; Q1=s.Q; R1=s.R; S1=s.S; fs=s.fs; nn=s.nn;
if nn(2)==0
   error('lss.balreal requires a system with at least one input')
end
P=gram(s,'c');
Q=gram(s,'o');
[Up,Sp,Vp]=svd(P);
Rp=(Up*diag(sqrt(diag(Sp))))';
[Uq,Sq,Vq]=svd(Q);
Rq=(Uq*diag(sqrt(diag(Sq))))';
[u, s1, v] = svd (Rq*Rp');
s1 = diag (s1);
ind=find(s1~=0);
sinv=zeros(size(s1));
sinv(ind)=1./s1(ind);
ssemi = diag (sqrt(sinv));
T = ssemi*u'*Rq;
Ti = Rp'*v*ssemi;
S=T*P*T';
Ab=T*A*Ti; Bb=T*B; Cb=C*Ti;
if nn(3)>0 %Stochastic system
    sys=lss(Ab,Bb,Cb,D,T*Q1*T',R1,T*S1,fs);
else
    sys=lss(Ab,Bb,Cb,D,fs);
end
sys=inherit(sys,s,'Balanced realization');
sys.sysMC=mceval1('balreal',sys);
end

function S=ctrb(s)
%CTRB computes the controllability matrix
%    S=ctrb(s)

A=s.A; B=s.B;  nx=s.nn(1);
S=B;
AB=B;
for i=2:nx
   AB=A*AB;
   S=[S AB];
end
end

function O=obsv(s)
%OBSV computes the observability matrix
%    O=obsv(s)

A=s.A; C=s.C;  nx=s.nn(1);
O=C;
CA=C;
for i=2:nx
   CA=CA*A;
   O=[O;CA];
end
end

function P=gram(s,opt,method)
%GRAM computes the controllability (observability) Gramian
%   P=gram(s,'c')  returns the controllability Gramian (default)
%   P=gram(s,'o')  returns the observability Gramian
%
%   The controllability Gramian P is the solution to
%      P=APA'+BB'     for discrete time systems
%      AP+PA'+BB'=0   for continuous time systems
%   The observability Gramian is implicitely defined by Q=gram(A',C');
%
%   Example:
%     m=rand(lss([2 1 0 1]))
%     gram(m,'o')
%     gram(m,'c')
%     mb=balreal(m)
%     gram(mb,'o')
%     gram(mb,'c')
%
%   method=1: direct solution of the system of linear equations
%   method=2: iterative solution
%   Default, method 1 is tried, and if it fails method two is applied

%$ Revision: 28-Oct-2019 $

if nargin<3; method=0; end
if nargin<2; opt='c'; end

if method==0
   try
      P=gram(s,opt,1);
   catch
      P=gram(s,opt,2);
   end
   return
elseif method~=1 & method~=2
   error('LSS.GRAM: method must be 0,1 or 2')
end

A=s.A; B=s.B; C=s.C;  nx=s.nn(1);
if strncmpi(opt,'o',1)
     B=C'; A=A';
elseif strncmpi(opt,'c',1)
     %
else
    error('Unknown option, ''c'' or ''o'' are allowed as second argument')
end
if s.fs>0 %discrete time
    if method==2
    % Stationary P
        maxiter=1000;
        P=eye(size(A));
        Pold=zeros(size(A));
        n=0;
        while norm(P-Pold)>1e-6 & n<maxiter
            n=n+1;
            Pold=P;
            P=A*P*A'+B*B';
            P=0.5*(P+P');
        end
        if n==maxiter; disp('Warning: maximum number of iterations attained in GRAM'), end
    else
        b=B*B';
        a=kron(A,A); % Use APA'=kron(A,A)vec(P)
        p=(eye(size(a))-a) \ b(:);
        P1=reshape(p,nx,nx);
        P=0.5*(P1+P1');
    end
else %continuous time system
    % Use vec(AXB)=kron(B',A) vec(X)
    b=B*B';
    a1=kron(eye(nx),A);
    a2=kron(A,eye(nx));
    p=-(a1+a2) \ b(:);
    P1=reshape(p,nx,nx);
    P=0.5*(P1+P1');
    A*P+P*A'+B*B';
end
if any(any(isnan(P))) | any(any(isinf(P)))
   error('LSS.GRAM: No solution found, system is probably unstable')
end
end

function [K,P,e]=lqe(s,varargin)
%LQE solves the continuous time stationary Riccati equation
%   [K,P,e]=lqe(s,Property1,Value1,...)
%   The stationary Ricatti equation is defined by
%   AP + PA' - PC'R^(-1)CP + Q = 0
%   K = PC'R^(-1)
%   where Q is factorized as Bv Qbar Bv'
%   A,B,C are the state space matrices of the model,
%   and Q, R, S are the covariance matrices of the noise sources.
%
%   K is the stationary Kalman gain
%   P is the stationary covariance matrix for the estimation error
%   e is the remaining error e=norm(AP + PA' - PC'R^(-1)CP + Q)
%
%   Property  Value         Description
%   -----------------------------------------------------------------
%   method    {0}|1|2|3     1 Eigenvector approach
%                           2 Iterative solution by solving ODE for P
%                           3 Algebraic solution when C is invertible
%                           0 Try 1 first, and catch with 2
%   disp      {'off'}|'on'  Display iteration info for method 2
%   maxiter   {1000}        Maximum number of iterations for method 2
%

opt=struct('method',0,'disp','off','maxiter',500);
opt=optset(opt,varargin);

if opt.method==0
   [K,P,e]=lqe(s,varargin{:},'method',1);
   if e>1e-8
      [K,P,e]=lqe(s,varargin{:},'method',2);
   end
   return
elseif opt.method~=1 & opt.method~=2
   error('LSS.LQE: method must be 0,1, 2 or 3')
end

A=s.A; B=s.B; C=s.C;
Q=s.Q; S=s.S; R=s.R;
nx=s.nn(1); ny=s.nn(4);
fs=s.fs;
if any(any(S~=0))
    error('LSS.LQE is not working for S~=0')
end
if all(all(Q==0)) & all(all(R==0))
    error('LSS.LQE works only for stochastic systems (Q~=0, R~=0)')
end

if fs>0
    error('LSS.LQE does not apply for discrete time systems, use DLQE instead')
end
if rank(R)<ny
    error('R must have full rank')
end

if opt.method==1
   H=[A' -C'*inv(R)*C;-Q -A];
   [V,D]=eig(H);
   V11=V(1:nx,1:nx);
   V21=V(nx+1:end,1:nx);
   P=V21*inv(V11);
   P=real(P);
   K=P*C'*inv(R);
   e=norm(A*P + P*A' - P*C'*inv(R)*C*P + Q);
elseif opt.method==2
   [T,tau]=timeconstant(A,NaN);
   P=zeros(nx);
   Pn=Q;
   maxiter=opt.maxiter;
   iter=0;
   M=C'*inv(R)*C;
   e=1;
   if strcmp(opt.disp,'on')
      disp('      Iteration   error  bisections iter bisection error')
   end
   while norm(e)>1e-10*norm(Pn) & iter<maxiter
      iter=iter+1;
      P=Pn;
      e=A*P + P*A' - P*M*P + Q;
      if any(any(isnan(e))) | any(any(isinf(e)))
          error('LSS.LQE: Inf or NaN obtained in P iterations')
      end
      mu=tau/100;   % Initial step size
      Pn=P+mu* e;
      Pn=0.5*(Pn+Pn');
      en=A*Pn + Pn*A' - Pn*M*Pn + Q;
      k=0; Pnn=0;
      while 0 & norm(en)>norm(e) & k<20  % Do not try bisection of step size
         k=k+1;
         mu=mu/2;
         Pnn=P+mu* e;
         Pnn=0.5*(Pnn+Pnn');
         en=A*Pnn + Pnn*A' - Pnn*M*Pnn + Q;
         if strcmp(opt.disp,'on')
           disp([iter norm(e) k norm(en)])
         end
      end
      if strcmp(opt.disp,'on')
         disp([iter norm(e)])
      end
      if k>0 & k<20  % bisection successful
         Pn=Pnn;
      end
      if any(any(isnan(Pn))) | any(any(isinf(Pn)))
          error('LSS.LQE: Inf or NaN obtained in P iterations')
      end
   end
   P=0.5*(Pn+Pn');
   K=P*C'*inv(R);
   if  norm(e)>1e-8
      disp(['LSS.LQE warning: remaining error ||AP+PA''-PC''R^(-1)CP+Q||=',num2str(norm(e))])
   end
   if iter==maxiter
      %error('LSS.LQE: maximum number of iterations used')
   end
elseif opt.method==3
   % solution for C invertible
   %  PMP-AP - PA' = Q
   %  (P-AM^-1) M (P-S^-1 A')  = Q+AM^-1A'
   %  P=AM^-1+(Q+AM^-1A')^0.5 M^(-0.5)
   M=C'*inv(R)*C;
   Minv=inv(M);  % M is not invertible generally
   RHS=Q+A*Minv*A';
   P=A*Minv+ sqrtm(RHS)*sqrt(Minv);
   K=P*C'*inv(R);
else
   error('LSS.LQE: method must be 1,2 or 3')
end
if any(any(isnan(P))) | any(any(isinf(P)))
   error('System must be stable to compute Gramians')
end
end


function [K,Pp,Pf]=dlqe(s)
%DLQE solves the discrete time stationary Riccati equation
%   [K,Pp,Pf]=dlqe(s)
%   A very simple iterative algorithm is used, nothing fancy at all,
%   to solve the stationary Ricatti equation
%   P=APA'-APC'(CPC'+R)^(-1)CPA'+Q
%   A,B,C are the state space matrices of the model,
%   and Q, R are the covariance matrices of the noise sources.
%   K is the stationary Kalman gain
%   Pp is the stationary covariance matrix for the prediction error
%   Pf is the stationary covariance matrix for the filtering error
%
%   Example:
%     m=exlti('motion1D')
%     [K,Pf,Pp]=dlqe(m)

A=s.A; B=s.B; C=s.C;
Q=s.Q; S=s.S; R=s.R;
nx=s.nn(1); fs=s.fs;
if any(any(S~=0))
    error('LSS.DLQE is not working when S~=0')
end
if all(all(Q==0)) & all(all(R==0))
    error('LSS.DLQE works only for stochastic systems (Q~=0, R~=0)')
end
Pp=eye(nx);
Ppold=zeros(nx);
n=0;
maxiter=1000;
if isnan(fs)
    error('LSS.DLQE does not apply for continuous time systems')
end
while norm(Pp-Ppold)>1e-10 & n<maxiter
    n=n+1;
    Ppold=Pp;
    L=Pp*C';
    K=L*inv(C*L+R);
    Pf=Pp-K*L';
    Pp=A*Pf*A'+Q;
    Pp=0.5*(Pp+Pp');
    Pf=0.5*(Pf+Pf');
end
if n==maxiter;
    disp('Warning: maximum number of iterations attained in DLQE'),
end
if any(any(isnan(Pp))) | any(any(isinf(Pp)))
   error('System must be stable to compute Gramians')
end
end

% =================================
% ----------Kalman filter----------
% =================================

function [x,V]=kalman(m,z,varargin)

%KALMAN implements various variants of the Kalman filter
%
%   [x,V]=kalman(m,z,Property1,Value1,...)
%
%   m      LSS object with model A,B,C, D, Q, R
%   z      SIG object with measurements
%   x      SIG object with state estimates
%          xhat=x.x and signal estimate yhat=x.y
%   V      Normalized sum of squared innovations, which should be a sequence
%          of chi2(nx) variables when the model is correct
%
%   Property  Value  Description
%   ------------------------------------------
%   alg       {1},2,3,4  1 stationary KF
%                        2, time-varying KF
%                        3, square root filter
%                        4, fixed interval KF smoother Rauch-Tung-Striebel
%                        5, sliding window KF, delivering xhat(t|y(t-k+1:t))
%                           where k is the length of the sliding window
%   k          k>0 {0}   Prediction horizon:
%                        0 for filter (default)
%                        1 for one-step ahead predictor,
%                        generally k>0 gives xhat(t+k|t) and y(t+k|t)
%                        for alg=1,2
%                        In case alg=5, k=L is the size of the sliding window
%   P0         {[]}      Initial covariance matrix
%                        Scalar value scales identity matrix
%                        Empty matrix gives a large identity matrix
%   x0         {[]}      Initial state matrix
%                        Empty matrix gives a zero vector
%   Q          {[]}      Process noise covariance (overrides the value in m.Q)
%                        Scalar value scales m.Q
%   R          {[]}      Measurement noise covariance (overrides the value in m)
%                        Scalar value scales m.R
%
%   Example:
%      %1. SISO
%      m=rand(lss([2 1 1 1],1));
%      u=getsignal('prbs',100);
%      z=simulate(m,u);
%      x1=kalman(m,z);
%      x2=kalman(m,z,'alg',2);
%      plot(z,x1,x2)
%      %2. Time series
%      m=rand(lss([2 0 1 1],1));
%      z=simulate(m,100);
%      x=kalman(m,z);
%      plot(z,x)  % Compare measurement with filtered C*xhat(t|t) version
%      xplot(x)   % Plot filtered state estimate
%      %3. Target simulation and tracking using the Kalman filter
%      m=exlti('ca2D');
%      z=simulate(m,10);
%      xhat=kalman(m,z);
%      xplot2(z,xhat,'conf',90,[2 3]);




opt=struct('alg',1,'k',0,'P0',[],'x0',[],'Q',[],'R',[]);
opt=optset(opt,varargin);
if isnan(m.fs)
    error('lss.kalman: only discrete time models')
end
[nx,nu,nv,ny]=size(m);
[A,B,C,D,Q,R,S,fs]=getss(m);
[N,nyz,nuz]=size(z);
if nyz~=ny
   error('LSS.KALMAN: Number of outputs ny must be the same in model and signal')
end
if nuz~=nu
   error('LSS.KALMAN: Number of inputs nu must be the same in model and signal')
end

y=z.y;
u=z.u;
if isempty(opt.P0)
    P0=1e6*eye(nx);
elseif isnumericscalar(opt.P0)
    P0=opt.P0*eye(nx);
elseif isnumeric(opt.P0) & size(opt.P0,1)==nx & size(opt.P0,2)==nx & norm(opt.P0-opt.P0')<1e-10
    P0=opt.P0;
else
   error('LSS.KALMAN: P0 must be empty, a scalar, or a symmetric (nx,nx) matrix')
end

if isempty(opt.x0)
    x0=zeros(nx,1);
elseif isnumeric(opt.x0) & size(opt.x0,1)==nx & size(opt.x0,2)==1
    x0=opt.x0;
else
   error('LSS.KALMAN: x0 must be empty or a (nx,1) vector')
end

if isempty(opt.Q)
    % Use Q from m
elseif isnumericscalar(opt.Q)
    Q=opt.Q*Q;
elseif iscov(opt.Q)  & size(opt.Q,1)==nx & size(opt.Q,2)==nx
    Q=opt.Q;
else
   error('LSS.KALMAN: Q must be empty, a scalar, or a covariance (nx,nx) matrix')
end

if isempty(opt.R)
    % Use R from m
elseif isnumericscalar(opt.R)
    R=opt.R*R;
elseif iscov(opt.R) & size(opt.R,1)==ny & size(opt.R,2)==ny
    R=opt.R;
else
   error('LSS.KALMAN: R must be empty, a scalar, or a covariance (ny,ny) matrix')
end

if isnumericscalar(opt.k) & opt.k>=0
   k=opt.k;
else
   error('LSS.KALMAN: k must be a scalar k>=0')
end
if opt.alg<1 | opt.alg>5 | opt.alg~=round(opt.alg)
   error('LSS.KALMAN: alg must be an integer in [1,5]')
end

if isempty(u)
   u=zeros(N,nu);  % Correct empty dimension
end


xhat=x0;
P=P0;
Px=zeros(N,nx,nx);
Py=zeros(N,ny,ny);
V=zeros(N,1);

if opt.alg==1;  % Stationary KF
  m.Q=Q; m.R=R;
  [K,Pp,Pf]=dlqe(m);
  % xfhat(t|t)=(A-KCA)xhat(t-1|t-1)+K(y(t)-Du(t))+(I-KC)Bu(t-1)
  % xphat(t+1|t)=(A-AKC)xhat(t|t-1)+AK(y(t)-Du(t))+Bu(t)
  % k=1, one step prediction
  [yhat,xhat]=dlsim(A-A*K*C,[K B-A*K*D],C,[zeros(ny,ny) D],[y u;zeros(1,ny+nu)],xhat);
  % xhat(k|k-1) and yhat(k|k-1) are computed
  xhat=xhat.';  % Simplifies matrix computations
  epsi=yhat(1:N,:)-y-u*D.';
  epsin=epsi*inv(C*Pp*C'+R);
  V=norm(epsi'*epsin);
  if k==0 % filtering
     xhat(:,end)=[];
     xhat=(eye(nx)-K*C)*xhat+K*y.'-K*D*u.';
     P=Pf;
  else
     xhat(:,1)=[];
     P=Pp;
  end
  if k>1 % Predict k steps ahead -> xhat(t+k|t), y(t+k|t)
    u=[u;zeros(k,nu)];
    for kk=2:k
      P=A*P*A' + Q;
      for t=size(xhat,2):-1:1
        xhat(:,t+1)=A*xhat(:,t) + B*u(t,:).';
      end
    end
  end
  if nu>0
     yhat=C*xhat+D*u.';
  else
     yhat=C*xhat;
  end
  Pxtmp(1,:,:)=P;
  Px=repmat(Pxtmp,[size(xhat,2),1,1]);
  Pytmp(1,:,:)=C*P*C'+R;
  Py=repmat(Pytmp,[size(xhat,2),1,1]);
  x=sig(yhat.',fs,u,xhat.',Py,Px);
end

if opt.alg==2  % time-varying KF
  yf=zeros(ny,N);
  yp=zeros(ny,N);
  xf=zeros(nx,N);
  xp=zeros(nx,N);
  xp(:,1)=xhat;
  Pfx=zeros(N,nx,nx);
  Pfy=zeros(N,ny,ny);
  for t=1:N,
    L=P*C';
    Pe=R+C*L;
    den=inv(Pe);
    epsi=y(t,:).'-C*xhat-D*u(t,:).';
    Epsi(:,t)=epsi;
    epsin=sqrtm(den)*epsi;
    Epsin(:,t)=epsin;
    K=L*den;
    V(t)=norm(epsin'*epsin);

    % Measurement update
    xhat=xhat+K*epsi;
    P=P-K*L';

    xf(:,t)=xhat; % xhat(t|t)
    Pfx(t,:,:)=0.5*(P+P.');
    yf(:,t)=C*xhat+D*u(t,:).';   % yhat(t|t)
    Pyt=C*P*C'+R;
    Pfy(t,:,:)=0.5*(Pyt+Pyt.');

    % Time update
    xhat=A*xhat+B*u(t,:).';
    P=A*P*A'+Q;
    P=.5*(P+P');

    xp(:,t)=xhat; % xhat(t+1|t)
    Ppx(t,:,:)=0.5*(P+P.');
    yp(:,t)=C*xhat+D*u(t,:).';   % yhat(t+1|t)
    Pyt=C*P*C'+R;
    Ppy(t,:,:)=0.5*(Pyt+Pyt.');

    if k>1 & t>k-1  % k-step prediction  % xhat(t+k|t)
      xhattmp(:,t)=xhat;
      Ptmp(t,:,:)=P;
      xhatout(:,t)=xhattmp(:,t-k+1);
      Ppx(t,:,:)=Ptmp(t-k+1,:,:);
      for kk=2:k
        xhatout(:,t)=A*xhatout(:,t) + B*u(t-kk+1,:)';
        Ppx(t,:,:)=A*squeeze(Px(t,:,:))*A' + Q;
      end
      yp(:,t)=C*xhatout(:,t)+D*u(t,:).';
      Pyt=C*squeeze(Ppx(t,:,:))*C'+R;
      Ppy(t,:,:)=0.5*(Pyt+Pyt.');
    end
  end
  if k==0
      x=sig(yf.',fs,u,xf.',Pfy,Pfx);
  else
      x=sig(yp.',fs,u,xp.',Ppy,Ppx);
  end
end

if opt.alg==3
  if k>1,
       error('lss.kalman: square root algorithm must have k=0 or k=1'),
  end
  yhat=zeros(ny,N);
  Epsi=zeros(ny,N);
  Epsin=zeros(ny,N);
  Ph=sqrtm(P);
  [U,S,Vtmp]=svd(Q);
  nv=rank(S);
  Bv=U(:,1:nv);
  Qh=sqrt(S(1:nv,1:nv));
  if k==1
    for t=1:N,
      qrmatrix = [sqrtm(R), C*Ph, zeros(ny,nv);...
                  zeros(nx,ny), A*Ph, Bv*Qh];
      [Qfact,Rfact]=qr(qrmatrix');
      Rfact=Rfact';
      Ph=Rfact(ny+1:ny+nx,ny+1:ny+nx);
      Sh=Rfact(1:ny,1:ny);
      Kh=Rfact(ny+1:ny+nx,1:ny);
      epsi=y(t,:).'-C*xhat -D*u(t,:).';
      epsin=inv(Sh)*epsi;
      xhat=A*xhat+Kh*epsin + B*u(t,:).';
      P=Ph*Ph';
      xhatout(:,t+1)=xhat;
      Px(t+1,:,:)=P;
      Py(t+1,:,:)=C*P*C'+R;
      yhat(:,t+1)=C*xhat  + B*u(t+1,:).';
      Epsi(:,t+1)=epsi;
      Epsin(:,t+1)=epsin;
      V(t+1)=norm(epsin'*epsin);
    end
    x=sig(yhat.',fs,u,xhatout.',Py,Px);
  elseif k==0
    A=eye(nx);
    B=zeros(nx,nu);
    for t=1:N,

      qrmatrix = [sqrtm(R), C*A*Ph, C*Bv*Qh;...
                  zeros(nx,ny), A*Ph, Bv*Qh];
      [Qfact,Rfact]=qr(qrmatrix');
      Rfact=Rfact';
      Ph=Rfact(ny+1:ny+nx,ny+1:ny+nx);
      Sh=Rfact(1:ny,1:ny);
      Kh=Rfact(ny+1:ny+nx,1:ny);

      % Include time update, xhat filtered estimate
      epsi=y(t,:).'-C*xhat -D*u(t,:)';

      epsin=inv(Sh)*epsi;
      xhat=xhat+Kh*epsin; % xhat(t|t)=xhat(t|t-1)+K*epsi
      P=Ph*Ph';
      xhatout(:,t)=xhat;
      Px(t,:,:)=P;
      Py(t,:,:)=C*P*C'+R;
      yhat(:,t)=C*xhat + D*u(t,:)';
      Epsi(:,t)=epsi;
      Epsin(:,t)=epsin;
      V(t)=norm(epsin'*epsin);
      % Time update
      xhat=A*xhat;     %  xhat(t+1|t)
      if ~isempty(u) & ~isempty(B);
        xhat=xhat + B*u(t,:)';
      end
    end
    x=sig(yhat.',fs,u,xhatout.',Py,Px);
  end
end

if opt.alg==4  % Smoothing
  if k~=0,
       disp('Warning: Kalman smoothing must have prediction horizon k=0'),
  end
  % Forward recursion
  xhatf=kalman(m,z,varargin{:},'alg',2,'k',0);
  xhatp=kalman(m,z,varargin{:},'alg',2,'k',1);
  xf=xhatf.x.';
  Pf=xhatf.Px;
  xp=xhatp.x.';
  Pp=xhatp.Px;
  xs=zeros(nx,N);
  Ps=zeros(N,nx,nx);
  ys=zeros(ny,N);
  Pys=zeros(N,ny,ny);
  xs(:,N)=xf(:,N);
  Ps(N,:,:)=Pf(N,:,:);
  % Backward recursion
  for t=N-1:-1:1;
    % x(t|N)=x(t|t)+P(t|t)A'(t)P(t+1|t) (x(t+1|N)-x(t+1|t))  RTS formulas
    xs(:,t)=xf(:,t) + squeeze(Pf(t,:,:))*A'*inv(squeeze(Pp(t+1,:,:)))*(xs(:,t+1)-xp(:,t+1));
    Ps(t,:,:)=squeeze(Pf(t,:,:)) + squeeze(Pf(t,:,:))*A'*inv(squeeze(Pp(t+1,:,:)))*(squeeze(Ps(t+1,:,:))-squeeze(Pp(t+1,:,:)))*inv(squeeze(Pp(t+1,:,:)))*A*squeeze(Pf(t,:,:));
    ys(:,t)=C*xs(:,t)+ D*u(t,:).';
    Pytmp=C*squeeze(Px(t,:,:))*C'+R;
    Pys(t,:,:)=Pytmp;
    epsin=Pytmp*(ys(:,t)-y(t,:).');
    V(t)=norm(epsin'*epsin);
  end
  x=sig(ys.',fs,u,xs.',Pys,Ps);
end

if opt.alg==5 % Sliding window
  yf=zeros(ny,N);
  xf=zeros(nx,N);
  Px=zeros(N,nx,nx);
  Py=zeros(N,ny,ny);
  for tc=k+1:N,
    xhat=x0;
    P=P0;
    for t=tc-k:tc
      L=P*C';
      S=R+C*L;
      den=inv(S);
      epsi=y(t,:).'-C*xhat -D*u(t,:)';
      K=L*den;
      % Measurement update
      xhat=xhat+K*epsi;
      P=P-K*L';
      % Time update
      if t<tc;
        xhat=A*xhat + B*u(t,:)';
        P=A*P*A'+Q;
        P=.5*(P+P');
      end
    end
    yf(:,tc)=C*xhat + D*u(tc,:)';
    xf(:,tc)=xhat; % xhat(t|t-k+1:t)
    Px(tc,:,:)=P;
    Py(tc,:,:)=C*P*C'+R;
  end
  x=sig(yf.',fs,u,xf.',Py,Px);
end
try
  x.xlabel=m.xlabel;
  x.ylabel=m.ylabel;
  x.ulabel=m.ulabel;
  x.name=m.name;
end
end

% =================================
% ----------LTI conversions--------
% =================================


function sys=lss2nl2(s)
%LSS2NLw converts linear state space model to non-linear model object
%   s=lss2nl2(s)
%
%   See also: nl.nl2lss, nl
%
%   Examples:
%   lss2nl(rand(lss([2 2 2 2])))
%   See also: nl2lss, ss2armax, lss2arx, ltf.ltf2lss


[nx,nu,nv,ny]=size(s);

f=[mat2str(s.A),'*x(1:',num2str(nx),',:)'];
h=[mat2str(s.C),'*x(1:',num2str(nx),',:)'];
if nu>0
  f=[f,'+',mat2str(s.B),'*u(1:',num2str(nu),',:)'];
  h=[h,'+',mat2str(s.D),'*u(1:',num2str(nu),',:)'];
end
h=[h,'+e(1:',num2str(ny),',:)'];
nv=rank(s.Q);
if nv>0
  [Uq,Dq]=svd(s.Q);
  Qbar=Dq(1:nv,1:nv);
  Bv=Uq(:,1:nv);
  f=[f,'+',mat2str(Bv),'*v(1:',num2str(nv),',:)'];
else
    Qbar=[];
end
sys=nl2(f,h,[nx nv nu ny 0],s.fs);
sys.x0=zeros(nx,1);
sys.pv=Qbar;
sys.pe=s.R;
sys.xlabel=s.xlabel;
sys.ylabel=s.ylabel;
sys.ulabel=s.ulabel;
sys.name=s.name;
end

function sys=lss2nl(s)
%LSS2NL converts linear state space model to non-linear model object
%   s=lss2nl(s)
%
%   See also: nl.nl2lss, nl
%
%   Examples:
%   lss2nl(rand(lss([2 2 2 2])))
%   See also: nl2lss, ss2armax, lss2arx, ltf.ltf2lss

[nx,nu,nv,ny]=size(s);
f=[mat2str(s.A),'*x(1:',num2str(nx),',:)'];
h=[mat2str(s.C),'*x(1:',num2str(nx),',:)'];
if nu>0
  f=[f,'+',mat2str(s.B),'*u(1:',num2str(nu),',:)'];
  h=[h,'+',mat2str(s.D),'*u(1:',num2str(nu),',:)'];
end
sys=nl(f,h,[nx nu ny 0],s.fs);
sys.x0=zeros(nx,1);
sys.pv=s.Q;
sys.pe=s.R;
sys.xlabel=s.xlabel;
sys.ylabel=s.ylabel;
sys.ulabel=s.ulabel;
sys.name=s.name;
end

function sys=lss2ltf(s)
%LSS2LTF converts state space model to transfer function
%   s=lss2ltf(s)
%
%   Converts A,B,C,D in state space model
%     x[k+1] = A x[k] + B u[k]
%       y[k] = C x[k] + D u[k]
%   to b,a in transfer function y[k] = b(q)/a(q) u[k]
%   In contrast to the function lss2ltf, this method works for MIMO systems
%
%   Any noise model is neglected in lss2ltf
%
%   Examples:
%   See also: ss2armax, lss2arx, ltf/ltf2lss

[nx,nu,nv,ny]=size(s);
if nv>0
%    error('LSS.LSS2LTF: stochastic systems (nv>0) cannot be converted to LTF')
end
if isempty(s.A) & isempty(s.D)
   sys=ltf([nx nx 0 nu ny]);
   return
end

[A,B,C,D,Q,R,S,fs]=getss(s);
if nx==0 % gain
   b(:,1,:)=D;
   a=1;
else
   a=real(poly(A));
   if nu>0
       for j=1:ny
           for i=1:nu
               b(j,:,i)=real(poly(A-B(:,i)*C(j,:))) + (D(j,i)-1)*a;
           end
       end
   else
      b=[];
   end
end
sys=ltf(b,a,s.fs);
sys=inherit(sys,s);
sys.sysMC=mceval1('lss2ltf',s);
end

function s=lss2armax(s)
%LSS2armax converts state space model to armax model
%   s=lss2armax(s)
%
%   Converts A,B,C,D,Q in state space model
%     x[k+1] = A x[k] + B u[k] + v[k]
%       y[k] = C x[k] + D u[k] + e[k]
%   to a,b,c in ARMAX model a(q) y[k] = b(q) u[k] +  c(q) e[k]
%
%   The first step is to convert to innovation form with one noise input
%     x[k+1] = A x[k] + B u[k] + Ke[k]
%       y[k] = C x[k] + D u[k] + e[k]
%   where K is the stationary Kalman gain computed by dlqe
%
%   If u=0, an ARMA model is computed.
%
%   Examples:
%   See also: lss2ltf, lss2arx, ltf/ltf2lss

% TBD

if isnan(s.fs)
    error('ss2armax only for discrete time models')
end
[nx,nu,nv,ny]=size(s);
[A,B,C,D,Q,R,S,fs]=getss(s);
a=real(poly(A));
if nu>0
    for j=1:ny
        for i=1:nu
            b(j,:,i)=real(poly(A-B(:,i)*C(j,:))) + (D(j,i)-1)*a;
        end
    end
else
    b=[];
end
if nv>0
    [G,D]=svd(s.Q);
    sigmav=sqrt(diag(D(1:nv)));
    sigmae=sqrt(diag(s.R));
    for j=1:ny
        for i=1:nv
            c(j,:,i)=real(poly(A-G(:,i)*C(j,:)));
        end
    end
else
    c=[];
end
sys=armax(a,b,c,s.fs);
sys=inherit(sys,s);
end

function Hf=lss2freq(s,varargin)
%LSS2FREQ converts the LSS model to a tranfer function frequency object
%   Hf=lss2freq(ltf,Property1,Value1,...)
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   N         {1024}      Number of frequency grid points
%   fmax      {'auto'}    Maximum frequency
%   f         {}          Frequency grid (overrides N and fmax)
%
%   Examples:
%     m=rand(lss([3 1 0 1]))
%     Hf=freq(m);      % Implicit call
%     Hf=lss2freq(m);   % Explicit call
%     bode(Hf);
%
%   See also:

opt=struct('N',1024,'f',[],'fmax','auto','MC',0);
opt=optset(opt,varargin);
fs=s.fs;
[nx,nu,nv,ny]=size(s);

if isempty(opt.f);
    if isnan(fs);  % Continuous time ltf
        if isnumericscalar(opt.fmax) & opt.fmax>0
            fmax=opt.fmax;
        else
            T=timeconstant(poly(s.A),fs);
            if isempty(T); T=10; end %ad-hoc
            fmax=8/T;
            % ad-hoc rule: Maximum interesting frequency is 8 times larger than
            % dominating pole and zero
        end
        f=linspace(0,1,opt.N)*fmax;
        f(1)=eps;
    else
        f=linspace(0,0.5,opt.N)*fs;
        f(1)=eps;
    end
else
    f=opt.f;
    f=f(:).';
end
N=length(f);
[A,B,C,D,Q,R,S]=getss(s);
nu=s.nn(2); ny=s.nn(4);
H=zeros(N,ny,nu);
if isnan(fs); % Continuous time
    for k=1:N
        H(k,:,:)=C*inv(i*2*pi*f(k)*eye(nx)-A)*B+D;
    end
else %discrete time
    for k=1:N
        H(k,:,:)=C*inv(exp(i*2*pi*f(k)/fs)*eye(nx)-A)*B+D;
    end
end
if ~isempty(s.sysMC)
    for k=1:length(s.sysMC)
        Hftmp=lss2freq(s.sysMC{k},varargin{:},'MC',0);
        HMC(k,:,:,:)=Hftmp.H;
    end
else
    HMC=[];
end
Hf=freq(H,f,fs,HMC);
Hf=inherit(Hf,s,'LSS -> FREQ');
end

function Phi=lss2spec(s,varargin)
%LSS2SPEC converts the stochastic part of an LSS model to a SPEC object
%   Hf=lss2spec(ltf,Property1,Value1,...)
%
%
%   Property   Value      Description
%   ---------------------------------------------------------
%   N         {1024}      Number of frequency grid points
%   fmax      {'auto'}    Maximum frequency
%   f         {}          Frequency grid (overrides N and fmax)
%
%   Examples:
%     m=rand(lss([3 0 1 1]))
%     Hf=spec(m);      % Implicit call
%     Hf=lss2spec(m);   % Explicit call
%     plot(Hf);
%
%   See also:

opt=struct('N',1024,'f',[],'fmax','auto','MC',0);
opt=optset(opt,varargin);
fs=s.fs;
[nx,nu,nv,ny]=size(s);
[A,B,C,D,Q,R,S]=getss(s);
if ~isempty(S) & any(S~=0)
   error('LSS.LSS2SPEC for S~=0 not yet implemented')
end
if nv>ny
   error('LSS.LSS2SPEC presumes nv<=ny')
end

if isempty(opt.f);
    if isnan(fs);  % Continuous time ltf
        if isnumericscalar(opt.fmax) & opt.fmax>0
            fmax=opt.fmax;
        else
            T=timeconstant(poly(s.A),fs);
            if isempty(T); T=10; end %ad-hoc
            fmax=8/T;
            % ad-hoc rule: Maximum interesting frequency is 8 times larger than
            % dominating pole and zero
        end
        f=linspace(0,1,opt.N)*fmax;
        f(1)=eps;
    else
        f=linspace(0,0.5,opt.N)*fs;
        f(1)=eps;
    end
else
    f=opt.f;
    if ~isnumeric(f), error('Frequency vector f must be numeric'), end
    f=f(:).';
end
N=length(f);
if nv==0
   if isempty(R) | all(all(R==0))
      error('LSS.LSS2SPEC is pointless when there is no noise (nv=0 and R=0)')
   else
      for k=1:N
         Phi(k,:,:)=R;
      end
   end
else
   [U,S,V]=svd(Q);
   B=U*sqrt(S);
   B=B(:,1:ny);
   if isnan(fs) % Continuous time
       for k=1:N
           H=C*inv(i*2*pi*f(k)*eye(nx)-A)*B;
           Phi(k,:,:)=H*H'+R;
       end
   else %discrete time
       for k=1:N
           H=C*inv(exp(i*2*pi*f(k)/fs)*eye(nx)-A)*B;
           Phi(k,:,:)=H*H'+R;
       end
   end
end
PhiMC=[];
if ~isempty(s.sysMC)
    for j=1:length(s.sysMC)
       Phitmp=lss2spec(s.sysMC{j},varargin{:},'MC',0);
       PhiMC(j,:,:,:)=Phitmp.Phi;
    end
else
   PhiMC=[];
end
Phi=spec(Phi,f,PhiMC);
Phi=inherit(Phi,s,'LSS -> SPEC');
end


function R=lss2covfun(s,varargin)
%LSS2COVFUN computes the LTI covariance function as a covariance object
%   R=lss2covfun(s,varargin)
%
%   Algorithm: Covariance function given by
%      R(0) = C * Pibar * C' + R
%      R(tau) = C * A^tau * Pibar * C' + C * A^(tau-1) * S,
%   where Pibar is the controllability Gramian (see gram).
%
%   Property   Value              Description
%   --------------------------------------------------------
%   taumax     {30}               Maximum lag for which the covariance
%                                 function is computed
%   MC         {length(s.sysMC)}  Number of Monte Carlo simulations to
%                                 compute confidence bound (0 means no bound)
%
%   Examples:
%     m=rand(lss([4 0 1 1]))
%     c=covfun(m);
%     plot(c)
%
%   See also:

opt=struct('MC',length(s.sysMC),'method','exact','taumax',30);
opt=optset(opt,varargin);

[nx,nu,nv,ny]=size(s);
[A,B,C,D,Q,R,S]=getss(s);
if nv==0
   if isempty(R) | all(all(R==0))
      error('LSS.LSS2COVFUN is pointless when there is no noise (nv=0 and R=0)')
   else
      r(1,:,:)=R;
      for k=2:opt.taumax
         r(k+1,:,:)=zeros(size(R));
      end
   end
else
   [U,D]=svd(Q);
   B=U*sqrt(D);
   P=gram(s,'o');
   r(1,:,:)=C*P*C' + R;
   r(2,:,:)=C*A*P*C'+C*S;
   for tau=2:opt.taumax
       r(tau+1,:,:)= C*A^tau*P*C'+C*A^(tau-1)*S;
   end
end
if ~isempty(s.sysMC)
    for i=1:length(s.sysMC)
       Rtmp=lss2covfun(s.sysMC{i},varargin{:},'MC',0);
       rMC(i,:,:,:)=Rtmp.R;
    end
else
   rMC=[];
end
R=covfun(r,0:opt.taumax,rMC);
R=inherit(R,s,'LSS -> COVFUN');
end


function r=rl(s,K)
%RL computes the root locus roots for feedback with -K of sys
%   r=rl(s,K)
%   K is a vector of gains
%   r is a (nx,length(K)) matrix with roots
[A,B,C,D,Q,R,S]=getss(s);
if size(D,1)~=size(D,2)
    error('Root locus only for square systems')
end
for j=1:length(K);
   rtmp=eig(A-B*inv(1/(eps+K(j))*eye(size(D))+D)*C);
   if j>1  % association of roots
       for i=1:length(rtmp);
           [dum,ind(i)]=min(abs(r(i,j-1)-rtmp));
       end
       r(:,j)=rtmp(ind);
   else
       r(:,j)=rtmp;
   end
end
end


function [z,p,k,zMC,pMC,kMC]=zpk(s)
%ZPK computes the zeros, poles and gain
%   [z,p,k,zMC,pMC,kMC]=zpk(s)
%   p is a row vector with poles
%   z is a (ny x nb x nu) matrix with zeros
%     for strictly proper transfer functions, z is filled up with Inf
%     according to the system theory definition G(Inf)=0 for proper systems
%   k is a (ny x nu) matrix with gains
%   zMC,pMC,kMC are the corresponding Monte Carlo arrays
%
%   Examples:
%     m1=rand(lss([3 1 0 1]))
%     [z,p,k]=zpk(m1)
%     spcodeoff
%     m2=rand(lss([2 2 0 2]))
%     [z,p,k]=zpk(m2)
%
%   See also: lss,ltf,ltf.zpk

% Using the LTF implementation.

[z,p,k,zMC,pMC,kMC]=zpk(ltf(s));
end



% =================================
% ----------LTI generators----------
% =================================

function sys=eye(s)
%EYE generates a unitary state space model of the same size as s.
%    Example: sys=eye(lss([1,1,0,2,2]));

[nx,nu,nv,ny]=size(s);
if nu~=ny
    error('EYE works only for square systems')
else
    sys=lss([],[],[],eye(nu),s.fs);
end
end


function sys=zeros(s)
%ZEROS generates a zero state space model of the same size as s.
%    Example: sys=zeros(lss([1,1,0,2,2]));

[nx,nu,nv,ny]=size(s);
sys=lss([],[],[],zeros(ny,nu),s.fs);
end

function sys=ones(s)
%ONES generates a all one state space model of the same size as s.
%    Example: sys=ones(lss([1,1,0,2,2]));

[nx,nu,nv,ny]=size(s);
sys=lss([],[],[],ones(ny,nu),s.fs);
end


function sys=rand(s,MC,varargin)
%RAND generates a random model or a cell array a random models
%   sys=rand(lss,MC,options)
%
%   sys is a cell array with MC samples if MC>1.
%   Default MC=1, in case sys is an LSS object (no cell array)
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
%   3. If s is a LSS model without uncertainty, then MC copies are returned
%
%   Example:
%      rand(lss([3,2,0,1]) generates a random model with 3 states,
%                        2 inputs, no noise and one output
%      rand(sys) generates a random model of the same size as sys
%                if sys is uncertain.
%
%   See also: ltf.rand, randpoly

if nargin<2; MC=1; end
if isempty(s.sysMC) & isempty(s.D) % Models from the prior
    nn=s.nn;
    fs=s.fs;
    nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
    for k=1:MC;
        [a,fs]=randpoly(nx,'fs',fs,varargin{:});
        if nu>0
            for i=1:nu;
                bi=randpoly(nx,'fs',fs,varargin{:});
                Gtmp=ltf(bi,a,fs);
                Gtmp=lss(Gtmp);
                A=Gtmp.A; B(:,i)=Gtmp.B; C=Gtmp.C; D(1,i)=Gtmp.D;
%                [A,B(:,i),C,D(1,i)]=ltf2lss(bi,a,'o');
            end
        else
            Gtmp=ltf([1 zeros(1,nx-1)],a,fs);
            Gtmp=lss(Gtmp);
            A=Gtmp.A; C=Gtmp.C;
%            [A,B,C,D]=ltf2lss(zeros(1,nx),a,'o');
            B=[]; D=[];
        end
        if nv>0
            for i=1:nv;
                gi=randpoly(nx-1,'fs',fs,varargin{:});
                Gtmp=ltf(gi,a,fs);
                Gtmp=lss(Gtmp);
                A=Gtmp.A; G(:,i)=Gtmp.B;
%                [A,G(:,i)]=ltf2lss(gi,a,'o');
            end
            Q=G*G';
            R=eye(ny);
        else
            Q=zeros(nx,nx); R=zeros(ny,ny); S=zeros(nx,ny);
        end
        if ny>1;
             %error('Only MISO systems implemented'),
             C=[C;rand(ny-1,nx)];
             D=[D;rand(ny-1,nu)];
        end
        sys{k} =lss(A,B,C,D,Q,R,fs);
    end
elseif length(s.sysMC)>0  % MC samples
    if MC<=length(s.sysMC);
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
end


% Stochastic systems

function su=uncertain(s,par,X,MC)
%UNCERTAIN introduces uncertainty in LSS coefficients par according to dist X
%   su=uncertain(s,par,X,MC)
%
%   par is a string containing the coefficients to randomize,
%       e.g. par='A(1,2), B(1,1)  C(2)' (the separators are arbitrary).
%   X   is a distribution of the same length (length(X)) as par
%       Use vertcat to specify multivariate distributions,
%       e.g. X=[udist(0,1); udist(0,1)]
%   MC  is the number of MC samples (default s.MC)
%   su  is the same as s, except for that the field sysMC has changed
%       and the mean values of the coefficients in par are set to E(X)
%
%   Example:
%     s=lss([-1 0;1 0],[1;0],[1 0]);  % DC motor
%     % 1. Uncertain time constant
%     dist1=udist(-1.5,-0.6);
%     su1=uncertain(s,'A(1,1)',dist1,10);
%     % 2. Uncertain time constant and gain
%     dist2=ndist([-1;1],0.1*eye(2));
%     su2=uncertain(s,'A(1,1) B(1)',dist2,10);  % Uncertain time constant
%     % Feedback with proportional controller
%     suc1=4*feedback(su1,3);
%     suc2=4*feedback(su2,3);
%     subplot(2,1,1), bodeamp(su,suc1,suc2)
%     subplot(2,1,2), plot(step(su),step(suc1),step(suc1))
%
%   See also: ltf.uncertain, lss.fix


if nargin<4; MC=s.MC; end
if ~isnumeric(MC) | length(MC)>1 | MC<0 | MC~=round(MC)
    error('LSS.UNCERTAIN: MC must be a positive integer')
elseif MC==0
   su=s;
   return
end
if isa(X,'pdfclass')
   ind=find(par=='(');
   ind2=find(par==')');
   N=length(ind);
   if N==length(X)
      su=s;                % Copy
      su.sysMC=rand(s,MC); % assure that suMC has correct dim
      su.MC=MC;
      x=rand(X,MC);        % random coefficients
      mu=E(X);
      for k=1:N
         cstr{k}=par(ind(k)-1:ind2(k));  % coefficient string
         eval(['su.',cstr{k},'=mu(k);'])
         for m=1:MC          % randomize the coefficients in sysMC
            eval(['su.sysMC{m}.',cstr{k},'=x(m,k);'])
         end
      end
   else
      error('LSS.UNCERTAIN: when setting a number of parameters to a PDFCLALSS, the length of X must match')
   end
else
   error('LSS.UNCERTAIN: X must be a PDFCLALSS object')
end
end


function s=fix(su)
%FIX removes uncertainty in LSS coefficients
%   s=fix(su)
%
%   This might be useful to speed up computations by avoiding
%   the Monte Carlo loop over uncertain parameters
%
%   See also: ltf.fix, lss.uncertain, lss.estimate

s=su;
s.sysMC=[];
end


function sys=estimate(s,sig,varargin)
%ESTIMATE estimates a state space model from the data in sig
%   sys=estimate(s,sig,varargin)
%
%   Uses the LTF method estimate.
%
%   Example:
%     G=rand(lss([2 1 0 1],1));
%     y=simulate(G,getsignal('prbs',100))+0.1*randn(100,1);
%     Ghat=estimate(G,y);
%     plot(G,Ghat)
%   See also: ltf.estimate, nl.estimate

sys=lss(estimate(ltf(s),sig,varargin{:}))
end

% =================================
% -------Data generators-----------
% =================================

function z=impulse(s,T)
%IMPULSE generates an impulse input and simulates the system
%   z=impulse(s,T)
%   This function computes analytically the impulse response at a time grid
%   For continuous time systems, the matrix exponential is used.
%   For discrete time systems, the Hankel coefficients are used
%
%   Argument    Description
%   --------------------------------------------------------
%   s           LSS object
%   T           Simulation time in seconds or samples.
%               Default it is  estimated from dominating pole
%   y           SIG object
%
%   Example:
%     G=lss(exlti('ltf2c'))
%     y=impulse(G);
%     plot(y)
%   See also: lss.step, lss.simulate

A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs; nn=s.nn;
nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
if nargin<2
    T=timeconstant(A,fs);
end
if fs>0 % discrete time
   if isnumericscalar(T)
      N=round(T*fs);
      if N>1e3; N=30; end  % Limit response for unstable systems
      if N<2
         error('LSS.IMPULSE: Simulation time T must be larger than 2/fs')
      end
      t=(0:N-1)'/fs;
   else
      error('LTF.IMPULSE: T must be either empty or a scalar T>0 for discrete systems')
   end
   if N<2
      error('LSS.IMPULSE: Simulation time T must be larger than 2/fs')
   end
   u=[ones(1,nu);zeros(N-1,nu)];
   t=(0:N-1)'/fs;
   [y,x]=dlsim(A,B,C,D,u);
else % continuous time system
   if isnumericscalar(T)
      N=1000;
      t=(0:N-1)'/N*T;
   elseif isvector(T) & all(T>=0)
      t=T(:);
      N=length(t);
   else
      error('LTF.IMPULSE: T must be either empty, a scalar T>0 or a time vector t>=0 for continuous time systems')
   end
    for k=1:length(t);
        E=real(expm([ [A B]*t(k); zeros(nu,nx+nu)]));
        Ak=E(1:nx,1:nx);
        x(:,k)=Ak*B*ones(nu,1);
        y(:,k)=C*x(:,k);
    end
    y=y.';
    x=x.';
    % Special continuous time representation
    if ~isempty(find(t==0))
       x=[zeros(2,nx);x];
       u=[zeros(1,nu);ones(1,nu);zeros(N,nu)];
       y=[zeros(2,ny);y];
       t=[t(1);t(1);t];
    else
       u=[zeros(N,nu)];  % Unit impulse
    end
end
% Recursive calls for Monte Carlo simulations
if ~isempty(s.sysMC)
    for i=1:length(s.sysMC);
        ztmp=impulse(s.sysMC{i},T);
        yMC(i,:,:)=ztmp.y;
        xMC(i,:,:)=ztmp.x;
    end
else
   yMC=[];
   xMC=[];
end
z=sig(y,t,u,x,yMC,xMC);
z.fs=fs;
z=inherit(z,s,'Impulse response');
if nargout==0
    if isinf(fs)
        plot(z,'conf',90)
    else
        staircase(z,'conf',90)
    end
end
end


function z=step(s,T)
%STEP generates a step input and simulates the system
%   y=step(s,T)
%
%   Argument    Description
%   --------------------------------------------------------
%   s           LSS object
%   T           Simulation time in seconds or samples,
%               or a time vector t for continuous time systems
%               Default T it is  estimated from dominating pole
%   y           SIG object
%
%   Example:
%     G=lss(exlti('ltf2c'))
%     y=step(G);
%     plot(y)
%   See also: lss.impulse, lss.simulate

A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs; nn=s.nn;
nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
if nargin<2
    T=timeconstant(A,fs);
end
if fs>0 % discrete time
   if isnumericscalar(T)
      if length(T)>1
         error('LSS.STEP: T cannot be a vector for discrete time systems')
      end
      N=round(T*fs);
      if N>1e3; N=30; end  % Limit response for unstable systems
      if N<2
         error('LSS.STEP: Simulation time T must be larger than 2/fs')
      end
      t=(0:N-1)'/fs;
   else
      error('LTF.STEP: T must be either empty or a scalar T>0 for discrete systems')
   end
   u=ones(N,nu);
   [y,x]=dlsim(A,B,C,D,u);
else
   if isnumericscalar(T)
      N=1000;
      t=(0:N-1)'/N*T;
   elseif isvector(T) & all(T>=0)
      t=T(:);
      N=length(t);
   else
      error('LTF.STEP: T must be either empty, a scalar T>0 or a time vector t>=0 for continuous time systems')
   end
    for k=1:length(t);
        E=real(expm([ [A B]*t(k); zeros(nu,nx+nu)]));
        Bk=E(1:nx,nx+1:nx+nu);
        x(:,k)=Bk*ones(nu,1);
        y(:,k)=C*x(:,k)+D*ones(nu,1);
    end
    y=y.';
    x=x.';
    % Special continuous time representation
    if ~isempty(find(t==0))
       x=[zeros(1,nx);x];
       u=[zeros(1,nu);ones(N,nu)];
       y=[zeros(1,ny);y];
       t=[t(1);t];
    else
       u=[ones(N,nu)];
    end
end

% Recursive calls for Monte Carlo simulations
if ~isempty(s.sysMC)
    for i=1:length(s.sysMC);
        ztmp=step(s.sysMC{i},T);
        yMC(i,:,:)=ztmp.y;
        xMC(i,:,:)=ztmp.x;
    end
else
   yMC=[];
   xMC=[];
end
z=sig(y,t,u,x,yMC,xMC);
z.fs=fs;
z=inherit(z,s,'Step response');
if nargout==0
    if isinf(fs)
        plot(z,'conf',90)
    else
        staircase(z,'conf',90)
    end
end
end


function k=dcgain(s)
%DCGAIN computes the DC gain of a system
if s.fs>0
   k=sum(s.b,2)./sum(s.a);
else
   if s.nn(3)==0
      k=s.b(:,1,:);
   elseif s.nn(3)>0
      k=zeros(s.nn(4),s.nn(5));
   else
      k=Inf*ones(s.nn(4),s.nn(5));
   end
end
end

function [z,xf]=simulate(s,in1,varargin)
%SIMULATE simulates a signal from an LSS model
%   [z,xf]=simulate(G,in1,Property1,Value1,...)
%
%
%   G is the LSS object
%   u is signal object with the same
%   sampling interval as the model (for instance, use u=sig(uvec,fs)).
%   If there are inconsistent sampling intervals, the one in the SIG object
%   has precedence.
%   xf is the final state, and can be forwarded as argument 'xi'
%   for simulation over different segments.
%
%   1. For discrete time systems, the low-level filter functions is applied
%      to each input-output MIMO channel.
%
%      For discrete time models, the input in1 is
%        - in1=T simulation length [s] of data to simulate for stochastic systems
%        - a SIG object with the same sampling interval as the model
%          (for instance, use in1=sig(uvec,fs)).
%
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
%      For continuous time models, in1 is
%         - a signal object for deterministic input-output models
%         - the simulations time T for stochastic systems
%      Note that continuous time stochastic input-output systems
%      are not supported
%
%   For non-uniformly sampled inputs, only continuous time models apply.
%
%   Property   Value/{Default}   Description
%   --------------------------------------------------------------------------
%   MC         {s.MC}            Number of Monte Carlo simulations
%                                Default value inherited from model or signal.
%   xi         zeros(nx,1)       Initial state as a nx=max([na nb-1]) vector
%   Ts         {timeconstant(s)/200} Sampling interval for fast sampling

%   For non-uniformly sampled inputs, only continuous models and
%   signal objects in in1 can be used.

%   in1   Description
%   --------------------------------------------------------
%   u     Input for input-output models as a SIG object u=sig(u,t)
%         For discrete time models, also a (N,nu) matrix is accepted
%   T     Simulation time for noise models (no input)
%
%
%
%   Example:
%     u=[getsignal('impulsetrain',120,40);getsignal('csquare',120,60)];
%     G=lss(exlti('ltf2c'));
%     y=simulate(G,u);
%     plot(y)
%
%     t=[0 0 10 10 50]; yv=[0 1 1 2 2]; y=sig(yv,t); g=ltf([1 2 1],[1 1 1])
%     simulate(g,y)
%

if nargin<2; in1=1000; end
[nx,nu,nv,ny]=size(s);
[A,B,C,D,Q,R,S,fs]=getss(s);
if isa(in1,'sig')
    MC=max([length(s.sysMC) in1.MC]);
elseif ~isempty(s.sysMC)  % Monte Carlo systems has precedence to s.MC
    MC=length(s.sysMC);
else
    MC=s.MC;
end

opt=struct('MC',MC,'xi',0,'Ts',0);
opt=optset(opt,varargin);

u=[];
if size(opt.xi,1)~=nx
    if opt.xi==0
        opt.xi=zeros(nx,1);
    else
        error('Initial state x0 must be either 0 or a row vector of appropriate size')
    end
end
if fs>0
    %-------------------------------
    %---discrete time model-------
    %-------------------------------
    if isa(in1,'sig')
        u=in1.y;
        [N,nudata]=size(u);
        if nudata~=nu
           error('LSS.SIMULATE: Number of inputs in u not compatible with model')
        end
        if isnan(fs)
           error('LSS.SIMULATE: discrete time models are simulate only with discrete time signals')
        end
        N=size(u,1);
    elseif isnumericscalar(in1)
        N=in1;
        u=randn(N,nu);  % white noise input
    elseif isa(in1,'double')
        u=in1;
        [N,nudata]=size(u);
        %if N<nudata; error('Input matrix u should have more rows than columns'), end
        if nudata~=nu
           error('Number of inputs in u not compatible with model')
        end
    else
        error('Inappropriate argument in1. Should be either input or number of data.')
    end
    t=(0:N-1)'/fs;

    if nv>0
        if isa(s.pv,'pdfclass')
           v=rand(s.pv,nx,N);
        else
           v=randn(nx,N);
        end
        if isa(s.pe,'pdfclass')
           e=rand(s.pe,ny,N);
        else
           e=randn(ny,N);
        end
        w=[v;e];
        if isempty(S); S=zeros(nx,ny); end
        [Utmp,Dtmp]=svd([Q S;S' R]);
        w=Utmp*sqrt(Dtmp)*w;
        v=w(1:nx,:);
        e=w(nx+1:nx+ny,:);
    else
        v=zeros(nx,N);
        e=zeros(ny,N);
    end
    if nv>0 & nu>0
        [y,x]=dlsim(A,[B eye(nx) zeros(nx,ny)],C,[D zeros(ny,nx) eye(ny)],[u v.' e.'],opt.xi);
    elseif nu==0 % No input
        [y,x]=dlsim(A,[eye(nx) zeros(nx,ny)],C,[zeros(ny,nx) eye(ny)],[v.' e.'],opt.xi);
    elseif nv==0 % No noise
        [y,x]=dlsim(A,B,C,D,u,opt.xi);
    else %??
        error('LSS.SIMULATE: Either u or v needed as input')
    end
else
  %-------------------------------
  %---continuous time model-------
  %-------------------------------
  if nu>0 % in1 gives input signal
    if nv>0
       error('LSS.SIMULATE: continuous time stochastic IO systems not implemented (both u and v terms in LSS model)')
    end
    if  ~isa(in1,'sig')
      error('LSS.SIMULATE: For input output systems, second input argument must be a SIG object')
    else
      if opt.Ts==0 % Compute default value
         tau=timeconstant(s.A,fs);  % Dominating pole determines time grid
         Ts=tau/200; % 200 times time constant of system fast enough sampling
         Ts=min([Ts (in1.t(end)-in1.t(1))/200]);
      else
         Ts=opt.Ts;
      end
      [N,nudata]=size(in1.y);
      if nudata~=nu
        error('LSS.SIMULATE: Number of inputs in SIG object not compatible with model')
      end
      u=in1;
      tu=u.t;
      indstep=find(diff(tu)==0);      % Steps or impulses in input here
      indimp=find(diff(indstep)==1);  % Impulses in input here
      indbreak=indstep;               % Segment boundaries
      indbreak(indimp+1)=[];          % Remove second break point for each impulse
      indimp=indstep(indimp);
      if isempty(indbreak)
         indbreak=[1;N];
      end
      if indbreak(1)>1
         indbreak=[1; indbreak(:)];
      end
      if indbreak(end)<N
         indbreak=[indbreak(:); N];
      end
      xi=opt.xi;                            % Initial state
      t=[];
      ymat=[];
      xmat=[];
      umat=[];
      for k=1:length(indbreak)-1;           % Simulate one interval at the time
          t1=tu(indbreak(k));
          t2=tu(indbreak(k+1));
          Nk=ceil((t2-t1)/Ts);              % Nr of samples in segment k
          if Nk==1; error('STEP: sampling interval Ts too small'); end
          Tsk=(t2-t1)/Nk;                   % Sampling interval in segment k
          tk=linspace(t1,t2,Nk);
   %      tk(1)=[];            % Avoid segment boundaries with steps and impulses
          uk=u(indbreak(k):indbreak(k+1),:);% Segment k in u
          ukint=interp(uk,tk,'method','hold','degree',1); % Linear interpolation
          ukmat=ukint.y;                    % sig to matrix
          ukint=ukmat;                      % Input to fast filter
          if any(find(indbreak(k)==indimp))   % Fix for impulses
             impamp=uk.y(2,:)-uk.y(1,:);      % Impulse area
             ukint=[impamp/Tsk; ukint];
                  % Add extra sampel to ukint with a short pulse
          elseif any(find(indbreak(k)==indstep)) % Fix for steps
             %none
          end
          sdk=c2d(s,1/Tsk,'ZOH');           % Fast sampling of system
          [A,B,C,D,Q,R,S]=getss(sdk);
          [ykmat,xkmat]=dlsim(A,B,C,D,ukint,xi);
          if  any(find(indbreak(k)==indimp))       % Fix for impulses
             if length(t)>1 % Exception for first segment
                t=[t t(end) tk];                           % Repeat
                umat=[umat;umat(end,:)+impamp; ukmat];  % Repeat
                ymat=[ymat;ymat(end,:); ykmat(2:end,:)]; % Repeat
                xmat=[xmat;xmat(end,:); xkmat(2:end,:)]; % Repeat
             else
                t=[0 0 tk];
                umat=[uk.y(1:2,:);ukmat];
                %ymat=[ones(2,1)*ykmat(1,:); ykmat(2:end,:)];
                ymat=[zeros(2,ny); ykmat(2:end,:)];
                xmat=[ones(2,1)*xkmat(1,:); xkmat(2:end,:)];
             end
          elseif any(find(indbreak(k)==indstep))  % Fix for steps
             if length(t)>1 % Exception for first segment
                t=[t tk];                         % Repeat once
                umat=[umat;ukmat];              % Repeat once
                ymat=[ymat;ykmat];              % Repeat once
                xmat=[xmat;xkmat];              % Repeat once
             else
                t=[tk(1) tk];
                umat=[uk.y(1,:);ukmat];
                ymat=[ykmat(1,:);ykmat];
                xmat=[xkmat(1,:);xkmat];
             end
          else
             umat=[umat;ukmat];                % Save u matrix
             t=[t tk(1:end)];                  % and t vector
             ymat=[ymat;ykmat];                % and y matrix
             xmat=[xmat;xkmat];                % and x matrix
          end
          xi=xkmat(end,:).';
      end
      y=ymat; u=umat; x=xmat;
    end
  else  % noise only model
        if isnumericscalar(in1)
            T=in1;
            if opt.Ts==0 % Compute default value
                % Dominating pole determines time grid
                tau=timeconstant(s.A,fs);
                Ts=tau/200; % 200 times time constant of system
                Ts=min([Ts T/200]);
            else
                Ts=opt.Ts;
            end
            N=round(T/Ts);
            Ts=T/N;  % Adjust to get perfect fit to time interval T
            md=c2d(s,1/Ts,'ZOH');  % Fast sampling of system
            z=simulate(md,N);      % Fast simulation
            y=z.y; u=z.u; x=z.x; t=z.t;
        else
            error('For continuous time noise models, second input argument should be simulation time')
        end
  end
end

xf=x(end,:).';
% Recursive calls for Monte Carlo simulations
yMC=[];
xMC=[];

if opt.MC>0
   if ~isempty(s.sysMC)
       sMC=rand(s,opt.MC);
       for i=1:opt.MC;
           ztmp=simulate(sMC{i},in1,varargin{:},'MC',0,'Ts',opt.Ts);
           yMC(i,:,:)=ztmp.y;
           xMC(i,:,:)=ztmp.x;
       end
   elseif isa(in1,'sig') & ~isempty(in1.yMC)
       uMC=rand(in1,opt.MC);
       for i=1:opt.MC;
           ztmp=simulate(s,uMC{i},varargin{:},'MC',0);
           yMC(i,:,:)=ztmp.y;
           xMC(i,:,:)=ztmp.x;
       end
   else
       for i=1:opt.MC;
           ztmp=simulate(s,in1,varargin{:},'MC',0);
           yMC(i,:,:)=ztmp.y;
           xMC(i,:,:)=ztmp.x;
       end
   end
end
z=sig(y,t,u,x,yMC,xMC);
z.fs=fs;
z.MC=size(yMC,1);  % XXX why needed
z=inherit(z,s,'Simulation');
if ~isempty(z.name)
   z.name=['Simulation of ',z.name];
end

if nargout==0
    plot(z,'conf',90)
end
end



% =================================
% -------Utility functions---------
% =================================

function [nx,nu,nv,ny]=size(s,dim)
%SIZE returns the sizes nn=[nx,nu,nv,ny]
%   [nx,nu,nv,ny]=size(s)
%   Special cases:
%   n=size(s) gives n=[ny nu]. Useful in for instance s+ones(size(s))
%   ny=size(s,1);
%   nu=size(s,2);

nn = s.nn;
nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
if nargin==2
    if dim==1; nx=ny; end
    if dim==2; nx=nu; end
end
if nargout==1,
    na=[ny nu];
elseif nargout==0,
   disp(['nx = ',num2str(nx)])
   disp(['nu = ',num2str(nu)])
   disp(['nv = ',num2str(nv)])
   disp(['ny = ',num2str(ny)])
end
end

function texcode=tex(s,varargin)
%TEX creates latex code for LSS objects
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
%   filename   {''}             Name of the .tex file (none for '')
%   decimals   {1}              Number of decimals
%   env        {'align*'}    Tex environment, '' means no env
%
%   Examples:
%      tex(rand(lss([2 1 1 1])))
%
%   See also: textable, texmatrix, ltf.tex

opt=struct('filename','','decimals',1,'env','align*','format','%11.2g');
opt=optset(opt,varargin);
if ~isempty(opt.env)
    texcode=sprintf(['\\begin{',opt.env,'}']);
else
    texcode='';
end



A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs; nn=s.nn;
nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
if ~isempty(opt.env)
    texcode=sprintf(['\\begin{',opt.env,'}']);
else
    texcode='';
end
if nx>0
    if fs>0
        texcode=strvcat(texcode,'x(t+1) &=');
    else
        texcode=strvcat(texcode,'\dot{x}(t) &=');
    end
    texcode=strvcat(texcode,texmatrix(A,'env','','decimals',opt.decimals));
    texcode=strvcat(texcode,'x(t)');
    if nu>0 & any(B~=0)
        texcode=strvcat(texcode,' + ');
        texcode=strvcat(texcode,texmatrix(B,'env','','decimals',opt.decimals));
        texcode=strvcat(texcode,'u(t)');
    end
    if nv>0 & any(Q~=0)
        texcode=strvcat(texcode,' + ');
        texcode=strvcat(texcode,'v(t)');
    end
    texcode=strvcat(texcode,'\\');
end
texcode=strvcat(texcode,'y(t) &=');
if nx>0
    texcode=strvcat(texcode,texmatrix(C,'env','','decimals',opt.decimals));
    texcode=strvcat(texcode,'x(t)');
end
if nu>0  & any(D~=0)
    texcode=strvcat(texcode,' + ');
    texcode=strvcat(texcode,texmatrix(D,'env','','decimals',opt.decimals));
    texcode=strvcat(texcode,'u(t)');
end
if nv>0
    texcode=strvcat(texcode,' + ');
    texcode=strvcat(texcode,'e(t)');
    texcode=strvcat(texcode,' \\');
    if any(any(Q~=0))
        texcode=strvcat(texcode,'Q &=');
        texcode=strvcat(texcode,texmatrix(Q,'env','','decimals',opt.decimals));
        texcode=strvcat(texcode,' \\');
    end
    if any(any(S~=0))
        texcode=strvcat(texcode,'S &=');
        texcode=strvcat(texcode,texmatrix(S,'env','','decimals',opt.decimals));
        texcode=strvcat(texcode,' \\');
    end
    if any(any(R~=0))
        texcode=strvcat(texcode,'R &=');
        texcode=strvcat(texcode,texmatrix(R,'env','','decimals',opt.decimals));
        texcode=strvcat(texcode,' \\');
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
end;
end


function disp(m)
%DISP gives a ascii formatted printout of the LSS object
%   disp(sys)  or  sys
format='%11.2g';
if isnan(m.fs)
    %disp('Continuous time state space model.')
else
    fsstr=num2str(m.fs);
end
% disp('State Space model')
%disp([' A=' mat2str(m.A,4)])
%disp([' B=' mat2str(m.B,4)])
%disp([' C=' mat2str(m.C,4)])
%disp([' D=' mat2str(m.D,4)])
A=m.A; B=m.B; C=m.C; D=m.D;
Q=m.Q; R=m.R; S=m.S;
fs=m.fs;
nn=m.nn;

nx=nn(1); nu=nn(2); nv=nn(3); ny=nn(4);
if isempty(A) & isempty(B) & isempty(C) & isempty(D)
    disp(['Unspecified LSS model with ',num2str(nx),' states, ',...
          num2str(nu),' inputs, ',...
          num2str(nv),' process noise dimensions, and ',...
          num2str(ny),' outputs. '])
    return
end
if isnan(fs)
    strxn1='d/dt x(t) = ';
    strxn=' x(t)';
    strun=' u(t)';
    stryn='y(t)';
    strvn='v(t)';
    stren='e(t)';
else
    strxn1='x[k+1] = ';
    strxn=' x[k]';
    strun=' u[k]';
    stryn='y[k]';
    strvn='v[k]';
    stren='e[k]';
    end
streq=' = ';
    strplus=' + ';
if nx>0
    if nx>1
        left=['/ ';ones(nx-2,1)*'| ';'\ '];
        right=[' \';ones(nx-2,1)*' |';' /'];
        up=floor(nx/2);
        down=floor((nx-1)/2);
        uleft=left; uright=right;
    else
       left=''; right='';
       down=0; up=0;
       if nu<2
           uleft=left; uright=right;
       else
           uleft='('; uright=')';
       end
    end
    strxleft=[char(32*ones(up,length(strxn1)));strxn1;char(32*ones(down,length(strxn1)))];
    strx=[char(32*ones(up,5));strxn;char(32*ones(down,5))];
    stru=[char(32*ones(up,5));strun;char(32*ones(down,5))];
    strv=[char(32*ones(up,4));strvn;char(32*ones(down,4))];
    stry=[char(32*ones(up,4));stryn;char(32*ones(down,4))];
    stre=[char(32*ones(up,4));stren;char(32*ones(down,4))];
    streq=[char(32*ones(up,3));streq;char(32*ones(down,3))];
    strplus=[char(32*ones(up,3));strplus;char(32*ones(down,3))];
    str1=[ strxleft, left, num2str(A,format),right, strx];
    if nu>0
        str2=[ strplus, uleft, num2str(B,format), uright, stru];
    else
        str2='';
    end
    if nv>0
        str3=[ strplus, strv];
    else
        str3='';
    end
    str=[str1 str2 str3];
    str=strsimplify(lss,str);
    if size(str,2)<=60
       disp(str)
    elseif size(str,2)<=120
       ind=blankind(lss,str(:,1:60));
       disp(str(:,1:ind(end)))
       disp('')
       disp(str(:,ind(end)+1:end))
    else
       if isnan(fs)
          disp([strxn1,' A x(t) + B u(t) + v(t)'])
       else
          disp([strxn1,' A x[k] + B u[k] + v[k]'])
       end
    end
    disp('')
end
if ny>1
    left=[' /';ones(ny-2,1)*'| ';' \'];
    right=['\ ';ones(ny-2,1)*' |';'/ '];
    up=floor(ny/2);
    down=floor((ny-1)/2);
else
    if nx>1
        left='('; right=')';
    else
        left=''; right='';
    end
    down=0; up=0;
end
strx=[char(32*ones(up,5));strxn;char(32*ones(down,5))];
stru=[char(32*ones(up,5));strun;char(32*ones(down,5))];
stry=[char(32*ones(up,4));stryn;char(32*ones(down,4))];
stre=[char(32*ones(up,4));stren;char(32*ones(down,4))];
streq=' = ';   streq=[char(32*ones(up,3));streq;char(32*ones(down,3))];
strplus=' + ';   strplus=[char(32*ones(up,3));strplus;char(32*ones(down,3))];
str1=[stry,streq];
if nx>0
    str2=[ left, num2str(C,format),right, strx];
else
   str2='';
end

if nu>0
    if ny==1 & nu>1
        uleft='('; uright=')';
    else
        uleft=left; uright=right;
    end
    str3=[ strplus, uleft, num2str(D,format),uright, stru];
else
    str3='';
end
if ~all(all(R==0))
    str4=[ strplus, stre];
else
    str4='';
end
str=[str1 str2 str3 str4];
str=strsimplify(lss,str);
    if size(str,2)<=60
       disp(str)
    elseif size(str,2)<=120
       ind=blankind(lss,str(:,1:60));
       disp(str(:,1:ind(end)))
       disp('')
       disp(str(:,ind(end)+1:end))
    else
       if isnan(fs)
          disp(['  y(t) = C x(t) + D u(t) + e(t)'])
       else
          disp(['  y[k] = C x[k] + D u[k] + e[k]'])
       end
    end
    disp('')
if nv>0
    disp('')
    if nx>1
        left=[' /';ones(nx-2,1)*'| ';' \'];
        right=['\ ';ones(nx-2,1)*' |';'/ '];
        up=floor(nx/2);
        down=floor((nx-1)/2);
    else
       left=''; right='';
       down=0; up=0;
    end
    strQ='Q = Cov(v) = ';
    strQ=[char(32*ones(up,size(strQ,2)));strQ;char(32*ones(down,size(strQ,2)))];
    if any(any(Q~=0))
    str=[strQ,left,num2str(Q,format),right];
    if size(str,2)<=60
       disp(str)
    elseif size(str,2)<=120
       ind=blankind(lss,str(:,1:60));
       disp(str(:,1:ind(end)))
       disp('')
       disp(str(:,ind(end)+1:end))
    else
       disp(['  Q(t) too large to diplay'])
    end
    disp('')
    end
    strS='S = Cov(v,e) = ';
    strS=[char(32*ones(up,15));strS;char(32*ones(down,15))];
    if any(any(S~=0))
        disp([strS,left,num2str(S,format),right])
    end
end
if ~(all(all(R==0)))
    if ny>1
        left=['/ ';ones(ny-2,1)*' |';'\ '];
        right=['\ ';ones(ny-2,1)*' |';'/ '];
        up=floor(ny/2);
        down=floor((ny-1)/2);
    else
       left=''; right='';
       down=0; up=0;
    end
    strR='R = Cov(e) = ';
    strR=[char(32*ones(up,size(strR,2)));strR;char(32*ones(down,size(strR,2)))];
    if any(any(R~=0))
        disp([strR,left,num2str(R,format),right])
    end
end
disp('')
end



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
end

% =================================
% -------Plots---------------------
% =================================



function [A,B,C,D,Q,R,S,fs]=getss(s);
%GETLSS unpacks LSS object
%    [A,B,C,D,Q,R,S,fs]=getss(s);
A=s.A; B=s.B; C=s.C; D=s.D; Q=s.Q; R=s.R; S=s.S; fs=s.fs;
end

function [s1,s2]=lsscheck(s1,s2)
if isa(s1,'ltf')
    s1=lss(s1);
elseif isa(s1,'double')
    s1=lss([],[],[],s1,s2.fs);
%    s1=lss(0,0,0,s1);
elseif ~isa(s1,'lss')
    error('s1 must be LSS, LTF or a scalar')
end
if isa(s2,'ltf')
    s2=lss(s2);
elseif isa(s2,'double')
    s2=lss([],[],[],s2,s1.fs);
%    s2=lss(0,0,0,s2);
elseif ~isa(s2,'lss')
    error('s2 must be LSS, LTF or a scalar')
end

if isnan(s1.fs) & isnan(s2.fs)   % Continuous time
    fs=NaN;
elseif isnan(s1.fs) | isnan(s2.fs)
    error('LSS: cannot mix continuous and discrete time models')
elseif s1.fs==s2.fs
    fs=s1.fs;
else
    error('LSS: Multi-rate LTI operations are not possible')
end
end



% =================================
% -------Monte Carlo functions-----
% =================================

function sysMC=mceval1(op,s,varargin)
%MCEVAL1 utility function to compute MC data
if isa(s,'lss') & ~isempty(s.sysMC)
    for i=1:length(s.sysMC)
        sysMC{i}=feval(op,s.sysMC{i},varargin{:});
    end
else
    sysMC=[];
end
end

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
end

% =================================
% -------Help functions---------
% =================================

function ind=blankind(dummy,str)
%BLANKIND used to find column indeces of blanks in string matrices
[n,m]=size(str);
ind=0;
for k=1:m
   if str(:,k)==char(32*ones(n,1))
      ind=[ind k];
   end
end
end


function str=strsimplify(dummy,str)
%STRSIMPLIFY used to remove consecutive blanks in string matrices
ind=blankind(lss,str);
ind1=find(diff(ind)>1);
ind(1)=[];
indblank=[];
for k=1:length(ind1)-1
    blanklength=ind(ind1(k+1)-1)-ind(ind1(k));
    if blanklength>1
       indblank=[indblank ind(ind1(k)):ind(ind1(k+1)-1)-2];
    end
end
str(:,indblank)=[];
end


function [Q,R,S,fs]=lssdatacheck(dummy,A,B,C,D,varargin)
% Check data of lss object
% Either call it as:
%   lssdatacheck(A,B,C,D,fs)
% or
%   lssdatacheck(A,B,C,D,Q,R,fs)
% or
%   lssdatacheck(A,B,C,D,Q,R,S,fs)

if ~isnumeric(A) | ~isnumeric(B) | ~isnumeric(C) | ~isnumeric(D)
    error('A, B, C, and D must be numeric')
end
if ~isempty(D)
   [ny,nu] = size(D);
   [nxB,nuB] = size(B);
else
   ny=size(C,1); nu=0;
   D=zeros(ny,nu);
end
if isempty(A),
    nx=0; nx2=0;
    if ~isempty(B)
       error('If A is empty, B must be empty')
    end
    if ~isempty(C)
       error('If A is empty, C must be empty')
    end
else
    [nx,nx2] = size(A);
end
[nyC,nxC] = size(C);
[nyD,nuD]=size(D);

if  nx~=nx2
    error('A must be square')
end
if nu>0 & nxB~=nx
    error('Number of rows in B must equal the dimension of A')
end
if nxC~=nx
    error('Number of columns in C must equal the dimension of A')
end
if nx>0 & nu>0 & nyD~=ny
    % exception for systems with A=B=C=[] and noise systems B=D=[]
    error('Number of rows in D and C must be the same')
end
if nu>0 & (nuD~=nu)  % exception for systems with A=B=C=[]
    error('Number of columns in B must be equal to the number of columns in D')
end
Q=zeros(nx,nx);
R=zeros(ny,ny);
S=zeros(nx,ny);

if length(varargin)==0
    fs = NaN;
elseif length(varargin)==1
    fs = varargin{1};
elseif length(varargin)==3
    Q = varargin{1};
    R = varargin{2};
    fs = varargin{3};
    if ~iscov(Q)
        error('Q must be a valid covariance (positive semidefinite symmetric) matrix')
    end
    if ~iscov(R)
        error('R must be a valid covariance (positive semidefinite symmetric) matrix')
    end
elseif length(varargin)==4
    Q = varargin{1};
    R = varargin{2};
    S = varargin{3};
    fs = varargin{4};
    if ~iscov(Q)
        error('Q must be a valid covariance (positive semidefinite symmetric) matrix')
    end
    if ~iscov(R)
        error('R must be a valid covariance (positive semidefinite symmetric) matrix')
    end
    if  ~isnumeric(S)
        error('S must be numeric')
    end
else
    error('Incorrect number of arguments')
end
if ~isequal(size(Q),[nx nx])
   error('LSS: Q must be of size (nx,nx)')
end
if ~isequal(size(R),[ny ny])
    error('LSS: R must be of size (ny,ny)')
end
if ~isequal(size(S),[nx ny])
    error('LSS: S must be of size (nx,ny)')
end

if ~isnan(fs) & (~isnumericscalar(fs) | fs<0)
    error('fs must be a non-negative number or NaN')
end
end

end % methods
end % lss
