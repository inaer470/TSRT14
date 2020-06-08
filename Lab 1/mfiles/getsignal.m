function z=getsignal(ex,N,opt1,opt2);

%GETSIGNAL generates standard signals as SIG objects
%   z=getsignal(ex,N,opt1,opt2);
%   Standard examples of SIG objects of different properties are returned
%   depending on the value of ex. Optional N defines either
%   1. the number of samples for discrete time filters  (default 1024), or
%   2. the simulation time T for continuous time filters (default 10)
%
%   Ex         Description
%   -----------------------------------------
%   ones       A unit signal [1 1...1] with ny=opt1 dimensions
%   zeros      A zero signal [0 0...0] with ny=opt1 dimensions
%   pulse      A single unit pulse [1 0...0]
%   step       A unit step [1...1]
%   ramp       A unit ramp with an initial zero and N/10 trailing ones
%   square     Square wave with period length opt1
%   sawtooth   Sawtooth wave with period length opt1
%   pulsetrain Pulse train with period length opt1
%   prbs       Pseudo-Random Binary Sequency with basic period length opt1
%              (default N/100) and transition probability opt2 (default 0.5)
%   cones      A continuous time unit signal [1 1] with t=[0 N] and ny=opt1
%   czeros     A continuous time zero signal [0 0] with t=[0 N] and ny=opt1
%   impulse    A continuous time unit impulse
%   cstep      A continuous time step function [0 1 1] with t=[0 0 N]
%   csquare    Continuous time square wave with period length opt1
%   impulsetrain Impulse train with period length opt1
%   cprbs      Continuous time PRBS
%
%   sinc       sin(pi*t)/(pi*t) with t=k*T where T=opt1
%   diric      The periodic sinc function
%              sin(N*pi*t)/(N*sin(pi*t)) with t=k*T
%              where T=opt1 and N=opt2
%   gausspulse sin(pi*t)*p(t;sigma) with t=k*T where p
%              is the Gaussian pdf, T=opt1 and sigma=opt2
%   chirp1     sin(pi*(t+a*t^2) with t=k*T where T=opt1 and a=opt2
%
%   sin1       One sinusoid in noise
%   sin2       Two sinusoids in noise
%   sin2n      Two sinusoids in LP noise
%
%   Examples:
%     s=getsignal('chirp1');
%     subplot(2,1,1),    plot(s)
%     subplot(2,1,2),    plot(tfd(s))
%
%   See also: exlti, exltv, dbsig

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if nargin<2; N=1024; end
if ~isstr(ex)
   error('GETSIGNAL: example type must be a string')
end
switch lower(ex)
    case 'pulse'
        y=[1;zeros(N-1,1)];
        z=sig(y,1);
        z.name='Pulse signal';
        z.desc=['Example getsignal(''pulse'',',num2str(N),')'];
    case 'step'
        y=[ones(N,1)];
        z=sig(y,1);
        z.name='Step signal';
        z.desc=['Example getsignal(''step'',',num2str(N),')'];
    case 'ones'
        if nargin>2; nu=opt1; else nu=1; end
        y=[ones(N,nu)];
        z=sig(y,1);
        z.name='Unit signal';
        z.desc=['Example getsignal(''unit'',',num2str(N),',',num2str(nu),')'];
    case 'zeros'
        if nargin>2; nu=opt1; else nu=1; end
        y=[zeros(N,nu)];
        z=sig(y,1);
        z.name='Zero signal';
        z.desc=['Example getsignal(''zeros'',',num2str(N),',',num2str(nu),')'];
    case 'ramp'
        N1=round(0.9*N);
        y=[(0:N1-1)'/N1;ones(N-N1,1)];
        z=sig(y,1);
        z.name='Ramp signal';
        z.desc=['Example getsignal(''ramp'',',num2str(N),')'];
    case 'pulsetrain'
        if nargin<3; opt1=ceil(N/10); end  % Period length
        y=[1;zeros(opt1-1,1)];
        y=kron(ones(ceil(N/opt1),1),y);
        y=y(1:N);
        z=sig(y,1);
        z.name='Pulse train signal';
        z.desc=['Example getsignal(''pulsetrain'',',num2str(N),',',num2str(opt1),')'];
    case 'square'
        if nargin<3; opt1=ceil(N/5); end  % Period length
        y=[zeros(round(opt1/2),1);ones(round(opt1/2),1)];
        y=kron(ones(ceil(N/opt1),1),y);
        y=y(1:N);
        z=sig(y,1);
        z.name='Square wave signal';
        z.desc=['Example getsignal(''square'',',num2str(N),',',num2str(opt1),')'];
    case 'sawtooth'
        if nargin<3; opt1=round(N/10); end
        y=(0:opt1-1)'/opt1;
        y=kron(ones(ceil(N/opt1),1),y);
        y=y(1:N);
        z=sig(y,1);
        z.name='Sawtooth signal';
        z.desc=['Example getsignal(''sawtooth'',',num2str(N),',',num2str(opt1),')'];
    case 'prbs'
        if nargin<4; opt2=0.5; end
        if nargin<3; opt1=round(N/50); end
        if opt1<1 | opt1>N, error('Period length 1<opt1<N'), end
        if opt2>1 | opt2<0, error('Switching probability 0<opt2<1'), end
        Nu=ceil(N/opt1);
        u=cumprod(sign(rand(Nu,1)-opt2));
        y=kron(u,ones(opt1,1));
        z=sig(y(1:N),1);
        z.name='Pseudo-random binary signal';
        z.desc=['Example getsignal(''prbs'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
% Continuous time signals
    case 'impulse'
        if nargin<2; T=10; else, T=N; end
        y=[0;1;0;0];
        t=[0;0;0;T];
        z=sig(y,t);
        z.name='Impulse signal';
        z.desc=['Example getsignal(''impulse'',',num2str(T),')'];
    case 'cones'
        if nargin<2; T=10; else, T=N; end
        if nargin>2; nu=opt1; else nu=1; end
        y=[ones(2,nu)];
        t=[0 T];
        z=sig(y,t);
        z.name='Unit signal';
        z.desc=['Example getsignal(''cones'',',num2str(T),',',num2str(nu),')'];
    case 'czeros'
        if nargin<2; T=10; else, T=N; end
        if nargin>2; nu=opt1; else nu=1; end
        y=[zeros(2,nu)];
        t=[0 T];
        z=sig(y,t);
        z.name='Zero signal';
        z.desc=['Example getsignal(''czeros'',',num2str(T),',',num2str(nu),')'];
    case 'cstep'
        if nargin<2; T=10; else, T=N; end
        y=[0;1;1];
        t=[0;0;T];
        z=sig(y,t);
        z.name='Step signal';
        z.desc=['Example getsignal(''cstep'',',num2str(T),')'];
    case 'impulsetrain'
        if nargin<2; T=10; else, T=N; end
        if nargin<3; opt1=ceil(T/10); end  % Period length
        n=floor(T/opt1);
        y=repmat([0;1;0],n,1);
        y=[y;0];
        t=kron((0:opt1:T)',ones(3,1));
        t(end-1:end)=[];
        z=sig(y,t);
        z.name='Impulse train signal';
        z.desc=['Example getsignal(''impulsetrain'',',num2str(T),',',num2str(opt1),')'];
    case 'csquare'
        if nargin<2; T=10; else, T=N; end
        if nargin<3; opt1=ceil(T/5); end  % Period length
        n=floor(T/opt1);
        y=repmat([0;0;1;1],n,1);
        t=kron((0:opt1/2:T)',ones(2,1));
        t(1)=[];
        t(end)=[];
        z=sig(y,t);
        z.name='Continuous square wave signal';
        z.desc=['Example getsignal(''csquare'',',num2str(T),',',num2str(opt1),')'];
    case 'cprbs'
        if nargin<2; T=10; else, T=N; end
        if nargin<4; opt2=0.5; end
        if nargin<3; opt1=T/50; end
        if opt1<=0 | opt1>=T, error('Period length 0<opt1<T'), end
        if opt2>1 | opt2<0, error('Switching probability 0<opt2<1'), end
        Nu=ceil(T/opt1);
        u=cumprod(sign(rand(Nu,1)-opt2));
        y=kron(u,ones(2,1));
        t=kron((0:opt1:T)',ones(2,1));
        t(1)=[];
        t(end)=[];
        yy=y(2:2:end);
        ind=find(diff(yy)==0);
        y([2*ind 2*ind+1])=[];
        t([2*ind 2*ind+1])=[];
        z=sig(y,t);
        z.name='Pseudo-random binary signal';
        z.desc=['Example getsignal(''cprbs'',',num2str(T),',',num2str(opt1),',',num2str(opt2),')'];

    case 'sinc'
        if nargin<3; opt1=50/N; end
        t=(-floor((N-1)/2):ceil((N-1)/2))'*opt1;
        ind=find(t~=0);
        y=ones(N,1);
        y(ind)=sin(pi*t(ind))./(pi*t(ind));
        z=sig(y,t);
        z.name='Sinc signal sin(pi*t)/(pi*t)';
        z.desc=['Example getsignal(''sinc'',',num2str(N),',',num2str(opt1),')'];
    case 'diric'
        if nargin<3; opt1=5/N; end
        if nargin<4; opt2=20; end
        t=opt1*(-floor((N-1)/2):ceil((N-1)/2))';
        ind=find(t~=round(t));
        ind1=find(t/2==round(t/2));
        ind2=find((t+1)/2==round((t+1)/2));
        y(ind1)=ones(size(ind1));
        y(ind2)=-ones(size(ind2));
        y(ind)=sin(pi*opt2*t(ind))./(opt2*sin(pi*t(ind)));
        z=sig(y',t);
        z.name='Dirichlet signal sin(N*pi*t)/(N*sin(pi*t))';
        z.desc=['Example getsignal(''diric'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
    case 'gausspulse'
        if nargin<3; opt1=50/N; end   %T
        if nargin<4; opt2=10; end      %sigma
        t=opt1*(-floor((N-1)/2):ceil((N-1)/2))';
        y=1/sqrt(2*pi*opt2)*exp(-t.^2/opt2^2).*sin(pi*t);
        z=sig(y,t);
        z.name='Gauss pulse signal';
        z.desc=['Example getsignal(''gausspulse'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
    case 'chirp'
        if nargin<3; opt1=25/N; end   %f
        if nargin<4; opt2=0.08; end      %df
        t=(0:N-1)'*opt1;
        y=sin(pi*(opt2*t.^2));
        z=sig(y,1/opt1);
        z.name='Chirp signal';
        z.desc=['Example getsignal(''chirp'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
    case 'chirp1'
        if nargin<3; opt1=50/N; end   %T
        if nargin<4; opt2=0.08; end      %sigma
        t=(0:N-1)'*opt1;
        y=sin(pi*(opt2*t.^2));
        z=sig(y,1/opt1);
        z.name='Chirp signal';
        z.desc=['Example getsignal(''chirp'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
    case 'chirp2'
        if nargin<3; opt1=50/N; end   %T
        if nargin<4; opt2=0.8; end      %sigma
        t=(0:N-1)'*opt1;
        y=sin(pi*(opt2*t.^2));
        z=sig(y,1/opt1);
        z.name='Chirp signal';
        z.desc=['Example getsignal(''chirp'',',num2str(N),',',num2str(opt1),',',num2str(opt2),')'];
    case 'sin1'
        if nargin<3; opt1=30; end
        t=(0:N-1)';
        fs=2;
        f=0.4;
        y0=sin(2*pi*f/fs*t);
        V=ndist(0,1);
        y0=sig(y0,fs);
        y0.MC=opt1;
        z=y0+V;
        z.name='Sinusoid signal with noise';
        z.desc=['Example getsignal(''sin1'',',num2str(N),',',num2str(opt1),')'];
    case 'sin2'
        if nargin<3; opt1=30; end
        t=(0:N-1)';
        fs=2;
        f1=0.4;
        f2=0.35;
        e=randn(size(t));
        y0=sin(2*pi*f1/fs*t)+sin(2*pi*f2/fs*t);
        V=ndist(0,1);
        y0=sig(y0,fs);
        y0.MC=opt1;
        z=y0+V;
        z.name='Sinusoid signal with noise';
        z.desc=['Example getsignal(''sin2'',',num2str(N),',',num2str(opt1),')'];
    case 'sin2n'
        if nargin<3; opt1=30; end
        t=(0:N-1)';
        fs=2;
        f1=0.4;
        f2=0.35;
        e=randn(size(t));
        H=getfilter(2,0.25);
        e=filter(H.b,H.a,e);
        y0=sin(2*pi*f1/fs*t)+sin(2*pi*f2/fs*t)+e;
        y=y0+e;
        for k=1:opt1 % MC
           yMC(k,:)=y0+filter(H.b,H.a,randn(size(t)));
        end
        z=sig(y,t,[],[],yMC);
        z.name='Sinusoid signal with colored noise';
        z.desc=['Example getsignal(''sin2n'',',num2str(N),',',num2str(opt1),')'];
    otherwise
        error([ex,' is not a valid example']);
end
if nargout==0
    info(z)
    plot(z)
end
