% Experiment setup
N=16;
fs=1; T=1/fs;                 % Sampling frequency/interval
m0=[0 1 1 1 1 0 1 -1 0 -1 -1 -1 -1 0 -1 1]';  % True sensor positions
p0=[0.9 -0.2]';               % Initial target position
M=length(m0)/2;               % Number of sensors

% Sensor model
smod=exsensor('toa',M);       % TOA model
smod.x0=p0;                   % Initial state
smod.th=m0(:);                % Sensor grid
R=1e-3*eye(M);                % Measurement noise
smod.pe=R;

% Motion model
mmod=exmotion('ctpva2d',1/fs);

% Total model
tmod=addsensor(mmod,smod);
tmod.x0=[p0' 2*pi/N  -0.01 pi/2 2*pi/N]'; % Initial state
pv=tmod.pv;
tmod.pv=[];  % No process noise in simulation

% Plots
y=simulate(tmod,2*N);         % Two laps
xplot2(y)
hold on
plot(smod,'linewidth',2,'fontsize',18)
axis([-1.2 1.2 -1.2 1.2])
