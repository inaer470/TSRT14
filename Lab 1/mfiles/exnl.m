function m=exnl(ex,opt1)
%EXNL generates NL objects
%   m=exnl(ex);
%
%   Ex          Description
%   ----------------------------------------------
%   ctcv2d      Coordinated turn model, cartesian velocity, 2D, Ts=opt1
%   ct, ctpv2d  Coordinated turn model, polar velocity,     2D, Ts=opt1
%   cv2d        Cartesian velocity linear model in 2D, Ts=opt1
%   pfex        Classic particle filter 1D example used in many papers
%   vdp         Van der Pol system
%   ball        Model of a bouncing ball
%   pendulum    One-dimensional continuous time pendulum model
%   pendulum2   Two-dimensional continuous time pendulum model
%   ipfr        Isothermal plug flow reactor
%   nipfr       Non-isothermal plug flow reactor
%   ottmar      Another chemical example
%   auv
%   sensornetwork  Three fixed sensors and one moving platform,
%               opt1 picks out a motion model ('ctcv2d' default)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $




if length(ex)>3 & strcmp(lower(ex(1:4)),'list');
m={'ctcv2d','ctpv2d','cv2d','pfex','vdp','ball','pendulum','pendulum2'};
   return

elseif length(ex)>3 & strcmp(lower(ex(1:4)),'test');
   mtype=exnl('list')
   for k=1:length(mtype)
      disp(mtype{k})
      mm=exnl(mtype{k});
   end
   return
end

switch lower(ex)

case {'imugps2d'}
if nargin<2
   T=0.5;
else
   T=opt1;
end
nn=[5 3 2 0];
f1=['x(1,:)+',num2str(T),'*x(3,:)'];
f2=['x(2,:)+',num2str(T),'*x(4,:)'];
f3=['x(3,:)+',num2str(T),'*u(1)*cos(x(5,:))'];
f4=['x(4,:)+',num2str(T),'*u(2)*sin(x(5,:))'];
f5=['x(5,:)+',num2str(T),'*u(3)'];
h1='x(1,:)';
h2='x(2,:)';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
h=['[',h1,';',h2,']'];
m=nl(f,h,nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
m.pv=B*diag([0.1 0.1 0.01])*B';
m.pe=diag([10 10]);
m.px0=diag([10 10 100 10 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v1','v2','w'};
m.ylabel={'x1','x2'};
m.name='2D inertial model with acceleration and yaw rate inputs and position measurement';

case {'cv','ctpv2d'}
if nargin<2
   T=0.5;
else
   T=opt1;
end
nn=[5 0 2 0];
f1=['x(1,:)+2*x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),'/2).*cos(x(4,:)+x(5,:)*',num2str(T),'/2)'];
f2=['x(2,:)+2*x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),'/2).*sin(x(4,:)+(eps+x(5,:))*',num2str(T),'/2)'];
f3='x(3,:)';
f4=['x(4,:)+x(5,:)*',num2str(T)];
f5='x(5,:)';
h1='sqrt(x(1,:).^2+x(2,:).^2)';
h2='atan2(x(2,:),x(1,:))';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
h=['[',h1,';',h2,']'];
m=nl(f,h,nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0;0 0;1 0;0 0;0 1];
m.pv=B*diag([0.1 0.01])*B';
m.pe=diag([1 0.003]);
m.px0=diag([10 10 100 10 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v','h','w'};
m.ylabel={'R','phi'};
m.name='Coordinated turn model with polar velocity';

case {'ctcv2d'}
if nargin<2
   T=0.5;
else
   T=opt1;
end
nn=[5 0 2 0];
f1=['x(1,:)+x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),')-x(4,:)./(eps+x(5,:)).*(1-cos(x(5,:)*',num2str(T),'))'];
f2=['x(2,:)+x(4,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),')+x(3,:)./(eps+x(5,:)).*(1-cos(x(5,:)*',num2str(T),'))'];
f3=['x(3,:).*cos(x(5,:)*',num2str(T),')-x(4,:).*sin(x(5,:)*',num2str(T),')'];
f4=['x(4,:).*cos(x(5,:)*',num2str(T),')+x(3,:).*sin(x(5,:)*',num2str(T),')'];
f5='x(5,:)';
h1='sqrt(x(1,:).^2+x(2,:).^2)';
h2='atan2(x(2,:),x(1,:))';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
h=['[',h1,';',h2,']'];
m=nl(f,h,nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
m.pv=B*diag([0.001 0.001 0.001])*B';
m.pe=diag([1 0.0001]);
m.fs=1/T;
m.xlabel={'x1','x2','v1','v2','w'};
m.ylabel={'R','phi'};
m.name='Coordinated turn model with Cartesian velocity';

case 'vdp'
f='[x(2,:);(1-x(1,:).^2).*x(2,:)-x(1,:)]';
h='x';
m=nl(f,h,[2 0 2 0]);
m.name='Van der Pol system';
m.x0=[2;0];

case 'vdpdisc'
% Ts=0.2 and Euler sampling
f='[x(1,:)+0.2*x(2,:);x(2,:)+0.2*((1-x(1,:).^2).*x(2,:)-x(1,:))]';
h='x';
m=nl(f,h,[2 0 2 0]);
m.fs=5;
m.name='Discretized van der Pol system (Euler Ts=0.2)';
m.x0=[2;0];
m.px0=10*eye(2);
m.pe=1*eye(2);
m.pv=0*eye(2);

case 'auv'
   f='[th(1) th(2) 0; th(3) th(4) 0; 0 1 0]*x+[th(5); th(6);0]*u';
%   f='[th(1) th(2); th(3) th(4)]*x+[th(5); th(6)]*u';
   h='x';
   nn=[3 1 3 6];
   m=nl(f,h,nn);
   m.pe=ndist([0;0;0],diag([0.001 0.01 1])); % corresponds to auv.ws
   m.name='AUV steering dynamics';
   m.th=[0.1273     1.1069    -0.4665    -2.6960     1.2888    -2.5942];
   m.th=[1.529 0.190 -44.970 -5.355 1.490 -35.586];
   m.xlabel={'v','r','psi'};
   m.ylabel=m.xlabel;
   m.ulabel='deltar';
   m.thlabel={'a11','a12','a21','a22','b1','b2'};

case 'cv2d'
   if nargin<2
      T=0.5;
   else
      T=opt1;
   end
   A=[1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
   B=[T^2/2 0; 0 T^2/2; T 0; 0 T];
   C=[1 0 0 0; 0 1 0 0];
   R=0.01*eye(2);
   nn=[4 0 2 0];
   m=nl('[1 0 0.5 0; 0 1 0 0.5; 0 0 1 0; 0 0 0 1]*x','[1 0 0 0; 0 1 0 0]*x',nn,1/T);
   m.pv=B*B';
   m.pe=R;
   m.xlabel={'X','Y','vX','vY'};
   m.ylabel={'X','Y'};
   m.name='Constant velocity motion model';

case 'cv2dbearing'
   if nargin<2
      T=0.5;
   else
      T=opt1;
   end
   A=[1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
   B=[T^2/2 0; 0 T^2/2; T 0; 0 T];
   C=[1 0 0 0; 0 1 0 0];
   R=0.001*eye(1);
   nn=[4 0 1 2];
   m=nl('[1 0 0.5 0; 0 1 0 0.5; 0 0 1 0; 0 0 0 1]*x','atan2(x(2,:)-th(2),x(1,:)-th(1))',nn,1/T);
   m.pv=0.01*B*B';
   m.pe=R;
   m.xlabel={'X','Y','vX','vY'};
   m.ylabel={'Bearing'};
   m.name='Bearings-only constant velocity motion model';

case 'cv2dradar'
   if nargin<2
      T=0.5;
   else
      T=opt1;
   end
   A=[1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
   B=[T^2/2 0; 0 T^2/2; T 0; 0 T];
   C=[1 0 0 0; 0 1 0 0];
   R=diag([0.0001 0.3]);
   nn=[4 0 2 2];
   h='[atan2(x(2,:)-th(2),x(1,:)-th(1));sqrt((x(1,:)-th(1)).^2+(x(2,:)-th(2)).^2)]';
   m=nl('[1 0 0.5 0; 0 1 0 0.5; 0 0 1 0; 0 0 0 1]*x',h,nn,1/T);
   m.x0=[0.8 0.8 0 0]';
   m.pv=B*B';
   m.pe=R;
   m.xlabel={'X','Y','vX','vY'};
   m.ylabel={'Bearing','Range'};
   m.name='Bearings-only constant velocity motion model';

case 'cp2dradar'
   if nargin<2
      T=1;
   else
      T=opt1;
   end
   nn=[2 0 2 0];
   m=nl('x','[sqrt(x(1,:).^2+x(2,:).^2);atan2(x(2,:),x(1,:))]',nn,1/T);
   m.x0=[0.8 0.8]';
   m.pv=1*eye(2);
   m.pe=diag([1 0.01]);
   m.xlabel={'X','Y'};
   m.ylabel={'Range','Bearing'};
   m.name='Random walk motion model with radar measurement';

case 'pfex'
   f='x(1,:)/2+25*x(1,:)./(1+x(1,:).^2)+8*cos(t)';
   h='x(1,:).^2/20';
   m=nl(f,h,[1 0 1 0]);
   m.fs=1;
   m.x0=5;
   m.px0=5;  % P0=cov(x0)
   m.pv=10;  % Q=cov(v)
   m.pe=1;   % R=cov(e)

case 'ball'
  f='[x(2,:);-th(1)*x(2,:).^2.*sign(x(2,:))-9.8.*sign(x(1,:))]';
  h='abs(x(1,:))';
  m=nl(f,h,[2 0 1 1]);
  m.x0=[1;0];
  m.th=1;
  m.name='Bouncing ball';
  m.xlabel={'Modified height','Modified speed'};
  m.ylabel='Height';
  m.thlabel='Air drag';

case 'ipfr'
   disp('  Example from rate1 on page 119 in')
   disp('     Introduction to chemical engineering computing')
   disp('     Bruce A. Finlayson')
   disp('     Wiley, 2006.')
   disp('  Example modified by removing nuisance concentrations CB and CC')
   disp('  Original example not identifiable of u and k, so linear term included')
   nn=[1 0 1 2];       % [nx nu ny nth]
   f='-(2*th(1)/th(2))*x - (2/th(2))*x.^2';
   h=inline('x','t','x','u','th');
   m=nl(f,h,nn);
   m.th=[0.1; 0.4];   % True parameters
   m.x0=2;               % True initial state
   m.pe=1e-10;            % Measurement noise covariance
   m.pe=[];            % Measurement noise covariance
   m.name='Isothermal plug flow reactor';
case 'nipfr'
   disp('  Example from rateSO2 on page 121 in')
   disp('     Introduction to chemical engineering computing')
   disp('     Bruce A. Finlayson')
   disp('     Wiley, 2006.')
   nn=[2 0 1 12];
   f='nipfr';
   h=inline('x(1)','t','x','u','th');
   m=nl(f,h,nn);
   m.th=[-50; 0.167; 2.2; -4.1; 673.2; 1.02e4;...
         -14.96; 11070; -1.331; 2331; -11.02; 11570];
   m.x0=[1; 673.2];
   m.pe=1e-10;            % Measurement noise covariance
   m.pe=[];            % Measurement noise covariance
   m.name='Non-isothermal plug flow reactor';

case 'reaction'
nn=[1 1 1 2];       % [nx nu ny nth]
f=inline('[-th(1)*exp(-th(2)/(8.31441*u(1,:)))*x(1,:); th(1)*exp(-th(2)/(8.31441*u(1,:)))*x(1,:)]','t','x','u','th');
f=inline('[-th(1)*10000*exp(-th(2)*10000/(8.31441*u(1,:)))*x(1,:)]','t','x','u','th');
h=inline('[x(1,:)]','t','x','u','th');
m=nl(f,h,nn);
m.th = [5 2]';
%m.x0=[1 0]';
m.x0=[1]';
m.pe=[];
m.ulabel='Temperature';
m.xlabel='Concentration';
m.ylabel='Measured concentration';
m.name='Temperature depending reaction';

case 'ottmar'
   nn=[3 0 3 3];       % [nx nu ny nth]
   f=inline('[-th(1)*x(1); th(1)*x(1)-th(2)*x(2); th(2)*x(2)-th(3)*x(3)]','t','x','u','th');
   h=inline('x','t','x','u','th');
   m=nl(f,h,nn);
   m.th = [0.00977527; 0.000503385; 0.00102519];
   m.x0=[10; 0; 0];
   m.pe=[];
   m.name='Ottmar';

case 'sensornetwork'
   if nargin<2; opt1='ctpv2d'; end
   mm=eval(['exnl(''',opt1,''');']);  % Motion model f
   h='[sqrt((x(1,:)-th(1)).^2+(x(2,:)-th(2)).^2);sqrt((x(1,:)-th(3)).^2+(x(2,:)-th(4)).^2);sqrt((x(1,:)-th(5)).^2+(x(2,:)-th(6)).^2)]'

   N=10;
   m=nl(mm.f,h,[mm.nn(1:2) 3 6]);
   m.x0=[0;1;1;0;0];     % Initial state
   m.th=[0 0 1 2 2 1]';  % True sensor positions
   m.fs=mm.fs;
   m.name='Sensor network';
   m.xlabel=mm.xlabel;
   m.ylabel={'Range 1','Range 2','Range 3'};
   m.thlabel={'pX(1)','pY(1)','pX(2)','pY(2)','pX(3)','pY(3)'}
case 'pendulum'
   f='[-9.8/th(1)*sin(x(2,:))-th(2)/th(1)^2*x(1,:);x(1,:)]';
   h='x(2,:)';
   m=nl(f,h,[2 0 1 2]);
   m.th=[1;0];
   m.x0=[0;0.2];
   m.fs=NaN;
   m.xlabel={'omega','theta'};
   m.ylabel={'theta'};
   m.thlabel={'l','b/m'};
   m.name='Continuous time pendulum in 1D';
case 'pendulum2'
   f1='x(3,:)';
   f2='x(4,:)';
   f3='sin(x(1,:)).*cos(x(1,:)).*x(4,:).^2 - 9.8/th(1)*sin(x(1,:))-th(2)*x(3,:)';
   f4='-2./tan(x(1,:)).*x(4,:).*x(3,:)';
   f=['[',f1,';',f2,';',f3,';',f4,']'];
   h='[x(1,:);x(2,:)]';
   m=nl(f,h,[4 0 2 2]);
   m.th=[1;0];
   m.x0=[0;0;0.1;0.2];
   m.fs=NaN;
   m.xlabel={'theta','phi','thetadot','phidot'};
   m.ylabel={'theta','phi'};
   m.thlabel={'l','b/m'};
   m.name='Continuous time pendulum in 2D';
case 'pendulum2x'
   f1='th(1)*x(3,:)';
   f2='th(1)*x(4,:).*sin(x(1,:))';
   f3='-th(1)*sin(x(1,:)).*cos(x(1,:)).*x(4,:).^2 - 9.8*sin(x(1,:))-th(2)*x(3,:)';
   f4='2*th(1)*cos(x(1,:)).*x(4,:).*x(3,:)';
   f=['[',f1,';',f2,';',f3,';',f4,']'];
   h='[x(1,:);x(2,:)]';
   g='[0 0;0 0;th(1) 0;0 th(1)*sin(x(1,:))]';   % Noise model!
   m=nl(f,h,[4 0 2 2]);
   m.th=[1;0.5];
   m.x0=[0.2;0;0;0.5];
   m.fs=NaN;
   m.xlabel={'theta','phi','thetadot','phidot'};
   m.ylabel={'theta','phi'};
   m.thlabel={'Length l','Damping b/m'};
   m.name='Continuous time pendulum in 2D';
   m.desc='Special case of eq on page 161 in Physics Handbook';


case 'dpendulum'
   if nargin<2
      T=0.5;
   else
      T=opt1;
   end
   f=['[x(1,:)+',num2str(T),'*(-9.8/th(1)*sin(x(2,:))-th(2)/th(1)^2*x(1,:));x(2,:)+',num2str(T),'*x(1,:)]'];
   h='x(2,:)';
   m=nl(f,h,[2 0 1 2]);
   m.th=[1;0.2];
   m.x0=[0;0.2];
   m.fs=1/T;
   m.xlabel={'omega','theta'};
   m.ylabel={'theta'};
   m.thlabel={'l','b/m'};
   m.name='Discrete time pendulum';
otherwise
   error(['Unknown option to exnl: ',ex])
end
