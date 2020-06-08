function m=exmotion(ex,opt)
%EXMOTION generates various motion models as NL objects without outputs
%   m=exmotion(ex,opt);
%
%   Sensor models from for instance exsensor.m can then be added to the
%   motion models with the method nl.addsensor
%
%   Ex           Description
%   ----------------------------------------------
%   ctcv2d       Coordinated turn model, cartesian velocity, 2D, Ts=opt
%   ct, ctpv2d   Coordinated turn model, polar velocity,     2D, Ts=opt
%   ctpva2d      Coordinated turn model, polar velocity, acc state, 2D, Ts=opt
%   cv2d         Cartesian velocity linear model in 2D, Ts=opt
%   imu2d        Dead-reckoning of acceleration and yaw rate, Ts=opt
%   imukin2d     Two-dimensional inertial model with aX, aY and wX as inputs
%   imukin2dbias As imukin2d but with 3 bias states for the inertial measurem.
%   imukin3d     Three-dimensional inertial model with a and w as the 6D input
%   imukin3dbias As imukin3d but with 6 bias states for the inertial measurem.

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

switch lower(ex)


case {'imu2d'}
if nargin<2
   T=0.5;
else
   T=opt;
end
nn=[5 3 0 0];
f1=['x(1,:)+',num2str(T),'*x(3,:)'];
f2=['x(2,:)+',num2str(T),'*x(4,:)'];
f3=['x(3,:)+',num2str(T),'*(cos(x(5,:))*u(1)-sin(x(5,:))*u(2))'];
f4=['x(4,:)+',num2str(T),'*(sin(x(5,:))*u(1)+cos(x(5,:))*u(2))'];
f5=['x(5,:)+',num2str(T),'*u(3)'];
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
m=nl(f,[],nn);
m.x0=[10 10 0 0 0]';
B=[0 0 0;  0 0 0; 1 0 0;  0 1 0;  0 0 1];
m.pv=B*diag([0.1 0.1 0.01])*B';
m.px0=diag([10 10 100 100 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v1','v2','Psi'};
m.ulabel={'aX', 'aY', 'wZ'};
m.name='2D inertial model with acceleration and yaw rate inputs and position measurement';

case {'cv','ctpv2d'}
if nargin<2
   T=0.5;
else
   T=opt;
end
nn=[5 0 0 0];
f1=['x(1,:)+2*x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),'/2).*cos(x(4,:)+x(5,:)*',num2str(T),'/2)'];
f2=['x(2,:)+2*x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),'/2).*sin(x(4,:)+(eps+x(5,:))*',num2str(T),'/2)'];
f3='x(3,:)';
f4=['x(4,:)+x(5,:)*',num2str(T)];
f5='x(5,:)';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
m=nl(f,[],nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0;0 0;1 0;0 0;0 1];
m.pv=B*diag([0.1 0.01])*B';
m.px0=diag([10 10 100 10 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v','h','w'};
m.name='Coordinated turn model with polar velocity';


case {'ctcv2d'}
if nargin<2
   T=0.5;
else
   T=opt;
end
nn=[5 0 0 0];
f1=['x(1,:)+x(3,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),')-x(4,:)./(eps+x(5,:)).*(1-cos(x(5,:)*',num2str(T),'))'];
f2=['x(2,:)+x(4,:)./(eps+x(5,:)).*sin((eps+x(5,:))*',num2str(T),')+x(3,:)./(eps+x(5,:)).*(1-cos(x(5,:)*',num2str(T),'))'];
f3=['x(3,:).*cos(x(5,:)*',num2str(T),')-x(4,:).*sin(x(5,:)*',num2str(T),')'];
f4=['x(4,:).*cos(x(5,:)*',num2str(T),')+x(3,:).*sin(x(5,:)*',num2str(T),')'];
f5='x(5,:)';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
m=nl(f,[],nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
m.pv=B*diag([0.001 0.001 0.001])*B';
m.px0=diag([10 10 50 50 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v1','v2','w'};
m.name='Coordinated turn model with Cartesian velocity';

case 'ctpva2d'
if nargin<2
   T=0.5;
else
   T=opt;
end
nn=[6 0 0 0];
f1=['x(1,:)+2*x(3,:)./(eps+x(6,:)).*sin((eps+x(6,:))*',num2str(T),'/2).*cos(x(5,:)+x(6,:)*',num2str(T),'/2)'];
f2=['x(2,:)+2*x(3,:)./(eps+x(6,:)).*sin((eps+x(6,:))*',num2str(T),'/2).*sin(x(5,:)+(eps+x(6,:))*',num2str(T),'/2)'];
f3=['x(3,:)+x(4,:)*',num2str(T)];
f4='x(4,:)';
f5=['x(5,:)+x(6,:)*',num2str(T)];
f6='x(6,:)';
f=['[',f1,';',f2,';',f3,';',f4,';',f5,';',f6,']'];
m=nl(f,[],nn);
m.x0=[10 10 5 0 0 0.1]';
B=[0 0;0 0;T^2/2 0;T 0;0 T^2/2;0 T];
m.pv=B*diag([0.1 0.01])*B';
m.px0=diag([10 10 100 10 10 1]);
m.fs=1/T;
m.xlabel={'x1','x2','v','a','h','w'};
m.name='Coordinated turn model with polar velocity and acceleration state';

case 'cp2d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(2);
   B=[T*eye(2)];
   nn=[2 0 0 0];
   h=[mat2strformat(A),'*x(1:2,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y'};
   m.name='2D constant position motion model';


case 'cv2d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(4);
   A(1:2,3:4)=T*eye(2);
   B=[T^2/2*eye(2);T*eye(2)];
   nn=[4 0 0 0];
   h=[mat2strformat(A),'*x(1:4,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','vX','vY'};
   m.name='2D constant velocity motion model';

case 'ca2d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(6);
   A(1:4,3:6)=T*eye(4);
   A(1:2,5:6)=T^2/2*eye(2);
   B=[T^3/3*eye(2);T^2/2*eye(2);T*eye(2)];
   nn=[6 0 0 0];
   h=[mat2strformat(A),'*x(1:6,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','vX','vY','aX','aY'};
   m.name='2D constant acceleration motion model';


case 'cj2d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(8);
   A(1:6,3:8)=T*eye(6);
   A(1:4,5:8)=T^2/2*eye(4);
   A(1:2,7:8)=T^3/3*eye(2);
   B=[T^4/4*eye(2);T^3/3*eye(2);T^2/2*eye(2);T*eye(2)];
   nn=[8 0 0 0];
   h=[mat2strformat(A),'*x(1:8,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','vX','vY','aX','aY','jX','jY'};
   m.name='2D constant jerk motion model';

case 'cv3d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(6);
   A(1:3,4:6)=T*eye(3);
   B=[T^2/2*eye(3);T*eye(3)];
   nn=[6 0 0 0];
   h=[mat2strformat(A),'*x(1:6,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','Z','vX','vY','vZ'};
   m.name='3D constant velocity motion model';

case 'ca3d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(9);
   A(1:6,4:9)=T*eye(6);
   A(1:3,4:6)=T^2/2*eye(3);
   B=[T^3/3*eye(3);T^2/2*eye(3);T*eye(3)];
   nn=[9 0 0 0];
   h=[mat2strformat(A),'*x(1:9,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','Z','vX','vY','vZ','aX','aY','aZ'};
   m.name='3D constant acceleration motion model';

case 'cj3d'
   if nargin<2
      T=0.5;
   else
      T=opt;
   end
   A=eye(12);
   A(1:9,4:12)=T*eye(9);
   A(1:6,7:12)=T^2/2*eye(6);
   A(1:3,10:12)=T^3/3*eye(3);
   B=[T^4/4*eye(3);T^3/3*eye(3);T^2/2*eye(3);T*eye(3)];
   nn=[12 0 0 0];
   h=[mat2strformat(A),'*x(1:12,:)'];
   m=nl(h,[],nn,1/T);
   m.pv=B*B';
   m.xlabel={'X','Y','Z','vX','vY','vZ','aX','aY','aZ','jX','jY','jZ'};
   m.name='3D constant jerk motion model';

case 'imukin2d'
nn=[5 3 0 1];
f1=['x(1,:)+th(1)*x(3,:)'];
f2=['x(2,:)+th(1)*x(4,:)'];
f3=['x(3,:)+th(1)*u(1)*cos(x(5,:))'];
f4=['x(4,:)+th(1)*u(2)*sin(x(5,:))'];
f5=['x(5,:)+th(1)*u(3)'];
f=['[',f1,';',f2,';',f3,';',f4,';',f5,']'];
m=nl(f,[],nn);
m.x0=[10 10 5 0 0.1]';
B=[0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
Q=diag([0.1 0.1 0.01]);
m.pv=B*Q*B';
m.px0=diag([10 10 100 10 1]);
m.th=1;
m.fs=1;
m.xlabel={'X','Y','vX','vY','psi'};
m.ulabel={'aX','aY','wX'};
m.thlabel={'Ts'};
m.name='2D inertial model';

case 'imukin2dbias'
nn=[8 3 0 1];
f1=['x(1,:)+th(1)*x(3,:)'];
f2=['x(2,:)+th(1)*x(4,:)'];
f3=['x(3,:)+th(1)*(u(1)+x(6,:)).*cos(x(5,:))'];
f4=['x(4,:)+th(1)*(u(2)+x(7,:)).*sin(x(5,:))'];
f5=['x(5,:)-th(1)*(u(3)+x(8,:))'];
f6=['x(6,:)'];
f7=['x(7,:)'];
f8=['x(8,:)'];

f=['[',f1,';',f2,';',f3,';',f4,';',f5,';',f6,';',f7,';',f8,']'];
m=nl(f,[],nn);
m.x0=[10 10 5 0 0.1 0 0 0]';
B=[zeros(2,6);eye(6)];
Q=diag([0.1 0.1 0.01 0*ones(1,3)]);
m.pv=B*Q*B';
m.px0=diag([10 10 100 10 1 0 0 0]);
m.fs=1;
m.th=1;
m.xlabel={'X','Y','vX','vY','psi','baX','baY','bwX'};
m.ulabel={'aX','aY','w'};
m.thlabel={'Ts'};
m.name='2D inertial model with bias states';

case 'imukin3d'
   f='imukin3d';  % Kinematic model
   Q= blkdiag(1e-2*eye(3), 1e-5*eye(3));  % Acc and gyro noise
   T=1;
   m=nl(f,[],[10 6 0 1]);
   m.fs=1/T;
   m.th=T;
%   rocketmodel.ylabel={'X','Y','Z'};
   m.ulabel={'aX','aY','aZ','wX','wY','wZ'};
   m.xlabel={'X','Y','Z','vX','vY','vZ','q1','q2','q3','q4'};
   m.thlabel={'Ts'};
   m.name='3D kinematic model';
case 'imukin3dbias'
   f='imukin3dbias';  % Kinematic model with bias
   Q= blkdiag(1e-2*eye(3), 1e-5*eye(3),zeros(6));  % Acc and gyro noise/offset
   T=1;
   m=nl(f,[],[16 6 0 1]);
   m.fs=1/T;
   m.th=T;
   m.ulabel={'aX','aY','aZ','wX','wY','wZ'};
   m.xlabel={'X','Y','Z','vX','vY','vZ','q1','q2','q3','q4','baX','baY','baZ','bwX','bwY','bwZ'};
   m.thlabel={'Ts'};
   m.name='3D kinematic model with bias states';
case 'list'
   m={'imu2d','ctpv2d','ctcv2d','ctpva2d','cp2d','cv2d','ca2d','cj2d','cv3d','ca3d','cj3d','imukin2d','imukin2dbias','imukin3d','imukin3dbias'};
case 'test'
   l=exmotion('list');
   for k=1:length(l);
      disp(l{k})
      m=exmotion(l{k});
   end
otherwise
   error(['Unknown option to exmotion: ',ex])
end
