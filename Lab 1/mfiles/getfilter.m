function G=getfilter(n,fc,varargin)

%GETFILTER approximates ideal transfer functions of type LP, HP, BS, BP
%   m=getfilter(n,fc,Property1,Value1,...)
%
%   The function computes both the continuous and discrete time filters.
%   If the sampling frequency fs=0, then only the continuous
%   time filter is returned.
%   If fs is provided, the cut-off frequencies must satisfy fc<fs/2
%   Input parameters
%   fc        cut-off frequency or vector of frequencies, to be normalized by sampling frequency
%   n         filter order
%
%   Property   Value       Description
%   ---------------------------------------------------------
%   type      {'LP'} | 'HP' | 'BP'
%   alg       {'butter'} | 'cheby1'
%   fs        {2}          Sampling frequency,
%                          fs=0 corresponds to a continuous filter.
%   ripple    {0.5}        Ripple in decibels for Chebyshev
%
%
%   Examples:
%     H1=getfilter(4,0.3,'alg','butter');
%     H2=getfilter(4,[0.3 0.7],'alg','cheby1','type','bp');
%     H3=getfilter(4,0.7,'type','hp');
%     lti.plot(H1,H2,H3,'plottype','plot'), grid
%
%   See also: tftool,

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

opt=struct('alg','butter','type','lp','ripple',0.5,'fs',2);
opt=optset(opt,varargin);
if strncmpi(opt.alg,'cheby',5)
    if strcmpi(opt.alg(end-1:end),'II') | strcmpi(opt.alg(end),'2')
        opt.alg='cheby2';
    else
        opt.alg='cheby1';
    end
end
if n<1 | round(n)~=n;
   error('GETFILTER: n must be a positive integer.'),
end
fs=opt.fs;

if fs>0 & (any(fc>fs/2) | any(fc<=0))
   error('GETFILTER: fc must be inside the interval [0 fs/2].'),
end

% step 1
if fs>0
% Pre-warping
    wc = 2*fs*tan(pi*fc/fs);
else
    wc=2*pi*fc;
end

% step 2
if strcmpi(opt.type,'bp') | strcmpi(opt.type,'bs')
    if length(wc)<2,
        warning('fc must include two frequencies for BP and BS'),
        wc=wc*[0.9,1.1];
    end
    w0=sqrt(wc(1)*wc(2));
    wc=wc(2)-wc(1);
else
    if length(wc)>1, error('fc must include only one frequency for LP and HP'), end
end
if strcmpi(opt.type,'hp') | strcmpi(opt.type,'bs')
    wc=1/wc;
end

% step 3
% Compute poles
if strncmpi(opt.alg,'butter',6)
    % Poles according to book 'Signalbehandling'
    pc=wc*exp(i*pi*(0.5+0.5/n:1/n:1.5-0.5/n));
    ac=real(poly(pc));
    bc(n+1)=ac(n+1);
elseif strncmpi(opt.alg,'cheby1',6)
    epsilon=opt.ripple;
    % Poles according to book 'Signalbehandling'
    alpha=1/epsilon+sqrt(1+1/epsilon^2);
    r1=0.5*wc*(alpha^(1/n)-alpha^(-1/n));
    r2=0.5*wc*(alpha^(1/n)+alpha^(-1/n));
    phi=pi*(0.5+0.5/n:1/n:1.5-0.5/n);
    pc=r1*cos(phi)+i*r2*sin(phi);
    % Poles according to book 'Digital Filters'
    a=sinh(1/n*asinh(1/epsilon));
    b=cosh(1/n*asinh(1/epsilon));
    pc=-wc*a*sin(pi/2/n*(2*(1:n)-1)) + i*wc*b*cos(pi/2/n*(2*(1:n)-1));
    ac=real(poly(pc));
    bc(n+1)=ac(n+1);
elseif strncmpi(opt.alg,'cheby2',6)
    % Poles and zeros according to book 'Digital Filters'
    epsilon=opt.ripple;
    ws=2*wc;  %?
    K=epsilon*cosh(n*acosh(ws/wc));
    zc=i*ws./cos(pi/2/n*(2*(1:n)-1));
    r=ws./sqrt(cosh(pi/2/n*(2*(1:n)-1))+cosh(1/n*asinh(K))^2-1);
    phi=1./tan(-cot(pi/2/n*(2*(1:n)-1))*coth(1/n*asinh(K)))+pi;
    pc=r.*exp(i*phi);
    ac=real(poly(pc));
    bc=real(poly(zc));
    bc=bc*ac(n+1)/bc(n+1);
zc,pc,r,phi
elseif strncmpi(opt.alg,'ellip',6)
    disp('Not yet implemented')
else
    error('Unknown filter algorithm')
end

% step 4
if strcmpi(opt.type,'hp') | strcmpi(opt.type,'bs') % Substitute s->1/s
    ac=ac(end:-1:1);
    bc=bc(end:-1:1);
    bc=bc/ac(1);
    ac=ac/ac(1);
end
if strcmpi(opt.type,'bp') | strcmpi(opt.type,'bs') % Substitute s->s^2+w0^2
    pp{1}=1;
    for k=1:n
        pp{k+1}=conv(pp{k},[1 0 w0^2]);
    end
    for k=1:n+1;
        pa(k,:)=ac(k)*[zeros(1,k-1) pp{n+2-k} zeros(1,k-1)];
    end
    ac=sum(pa);
%    bc=[zeros(1,n) bc(end:-1:1)];
    bc=[bc zeros(1,n)];
    bc=bc/ac(1);
    ac=ac/ac(1);
end
G=ltf(bc,ac); % Continuous time TF object

% step 5
if fs>0 % Discretize
    G=c2d(G,fs,'bil');
end
G.name=[opt.alg,' filter of type ',opt.type];

if nargout<1
    plot(G)
end
