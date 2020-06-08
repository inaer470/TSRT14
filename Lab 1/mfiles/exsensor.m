function s=exsensor(ex,M,N,nx,varargin)
%EXSENSOR is a database with sensor models for sensor networks and SLAM
%   s=exsensor(ex,M,N,nx)
%   M sensors and N targets (default 1)
%   nx denotes the state dimension of each target (default 2)
%   The rows in h are h((m-1)*N+n,:)=sensor(pt(n),ps(m)), m=1,..,M,  n=1,...,N
%   The position for target n pt(n) is assumed to be x((n-1)*nx+1:(n-1)*nx+2)
%   The position for sensor m ps(m) is assumed to be
%        ps(m)=th((m-1)*2+1:(m-1)*2+2)            for sensor networks
%        ps(m)=x(N*nx+(m-1)*2+1:N*nx+(n-1)*nx+2)  for SLAM
%   If slam is appended to the string in ex, then a SLAM model is obtained
%
%   Target locations s.x0 and sensor locations s.th (for SN) are
%   set to uniformly distributed in [0,1]x[0,1]. Change these if desired.
%
%
%   Options for ex:
%   'toa'    ||pt(n)-ps(m)||                 2D TOA as range measurement
%   'tdoa1'  ||pt(n)-ps(m)||+x(N*nx+1)       2D TDOA with bias state
%   'tdoa2'  ||pt(n)-ps(m)||-||pt(k)-ps(m)|| 2D TDOA as range differences
%   'doa'    atan2(pt(n),ps(m))              2D DOA as bearing measurement
%   'rss1'   th(n)+th(N+1)*10*log10(||pt(n,1:2)-ps(m)||)  RSS with parameters
%   'rss2'   x(n,N*nx+1)+x(n,N*nx+2)*10*log10(||pt(n)-ps(m)||)  RSS with states
%   'radar'                              2D combination of TOA and DOA above
%   'gps2d'                              2D position
%   'gps3d'                              3D position
%   'mag2d'                              2D magnetic field disturbance
%   'mag3d'                              3D magnetic field disturbance
%   'quat'   sum(q.^2)=1                 Quaternion constraint
%   '*slam'                              * is one of above, th as augmented states
%
%
%  l=exsensor('list') gives a cell array with all options

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2, M=1; end
if nargin<3, N=1; end
if nargin<4, nx=2; end
if ~isscalar(M), error('EXSENSOR: number of sensors must be a scalar'), end
if M~=round(M), error('EXSENSOR: number of sensors must be an integer'), end
if ~isscalar(N), error('EXSENSOR: number of targets must be a scalar'), end
if N~=round(N), error('EXSENSOR: number of targets must be an integer'), end
if ~isscalar(nx), error('EXSENSOR: nx must be a scalar'), end
if nx~=round(nx), error('EXSENSOR: nx must be an integer'), end


if length(ex)>4 & strcmp(lower(ex(end-3:end)),'slam');
     slam=1;
     postype='x';
     posx=',:';
else
     slam=0;
     postype='th';
     posx='';
end
ex=[ex,'   '];

if strcmp(lower(ex(1:3)),'toa');
  if slam
     name='TOA SLAM';
     posind=nx*N;
     nn=[nx*N+2*M 0 M*N 0];
     th=[];
  else
     name='TOA';
     posind=0;
     nn=[2*N 0 M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['sqrt((x(',num2str(nx*(n-1)+1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(nx*(n-1)+2),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2)']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
  %ylabel={'TOA'};

elseif strcmp(lower(ex(1:3)),'doa');
    if slam
     name='DOA SLAM';
     posind=nx*N;
     nn=[nx*N+2*M 0 M*N 0];
     th=[];
  else
     name='DOA';
     posind=0;
     nn=[nx*N 0 M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['atan2(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,'),x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,'))']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
  pe=1e-2*eye(nn(3));
%  ylabel={'Bearing'};

elseif strcmp(lower(ex(1:4)),'rss1');
  % RSS where nuisance parameters are in th
    if slam
     name='RSS1 SLAM with parametric model';
     posind=nx*N;
     posind2=nx*N+2*M;
     postype='x';
     nn=[2*N+2*M+M+1 0 M*N M+1];
     th=[];
  else
     name='RSS1 with parametric model';
     posind=0;
     posind2=2*M;
     postype='th';
     nn=[2*N 0 M*N 2*M+N+1];
     th=[ones(N,1);-2;rand(2*M,1)];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,[postype,'(',num2str(n+posind2),')+',postype,'(',num2str(N+1+posind2),')*5*log10(sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2))']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
  pe=1e-1*eye(nn(3));

elseif strcmp(lower(ex(1:4)),'rss2');
  % RSS where nuisance parameters are in x
    if slam
     name='RSS2 SLAM with parameters in state';
     posind=nx*N;
     posind2=nx*N+2*M;
     postype='x';
     nn=[N+1+nx*N+2*M 0 M*N 0];
     th=[];
  else
     name='RSS2 with parameters in state';
     posind=0;
     posind2=nx*N;
     postype='th';
     nn=[N+1+nx*N 0 M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['x(',num2str(n+posind2),')+x(',num2str(N+1+posind2),')*5*log10(sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2))']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];


elseif strcmp(lower(ex(1:4)),'maxw');
  % Maxwell where nuisance parameters are in x
    if slam
     name='Maxwell SLAM with parameters in state';
     posind=nx*N;
     posind2=nx*N+2*M;
     postype='x';
     nn=[N+1+nx*N+2*M 0 M*N 0];
     th=[];
  else
     name='Maxwell with parameters';
     posind=0;
     posind2=nx*N;
     postype='th';
     nn=[6 0 3*M 3*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['1/norm(x(1:3)-th(',num2str(m),'*3-3+[1 2 3]))^5*(3*x(1:3)*x(1:3)''-norm(x(1:3))^2*eye(4))*x(4:6)']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
h

elseif strcmp(lower(ex(1:4)),'list');
s={'toa','tdoa1','tdoa2','doa','radar','gps2d','gps3d','rss1','rss2',...
   'mag2d','mag3d','quat'};
   return

elseif strcmp(lower(ex(1:4)),'test');
   mtype=exsensor('list')
   for k=1:length(mtype)
      disp(mtype{k})
      for M=1:2
         s=exsensor(mtype{k},M);
      end
   end
   return

elseif strcmp(lower(ex(1:5)),'tdoa1');
    if slam
     name='TDOA1 SLAM  with range offset state';
     posind=nx*N+1;
     postype='x';
     nn=[nx*N+1+2*M 0 M*N 0];
     th=[];
  else
     name='TDOA1 with range offset state';
     posind=0;
     postype='th';
     nn=[nx*N+1 0 M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2) + x(',num2str(nx*N+1),',:)']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
  %ylabel={'TOA with offset'};

elseif strcmp(lower(ex(1:5)),'tdoa2');
    if slam
     name='TDOA2 SLAM with pairwise differences';
     posind=nx*N;
     postype='x';
     nn=[nx*N+2*M 0 M*N*(M-1)/2 0];
     th=[];
  else
     name='TDOA2 with pairwise differences';
     posind=0;
     postype='th';
     nn=[nx*N 0 M*N*(M-1)/2 2*M];
  end
  h=['['];
  R=zeros(nn(3));
  r=0; % row index
  for n=1:N
     for m=1:M
       for i=m+1:M
          r=r+1;
          h=[h,['sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2) - sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*i-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*i),posx,')).^2);']];
          R(r,m)=1; R(r,i)=-1;
       end
     end
  end
  h(end)=[];
  h=[h,']'];
  pe=1e-4*R*R';
  %ylabel={'Time difference'};

elseif strcmp(lower(ex(1:5)),'radar');
    if slam
     name='RADAR SLAM';
     posind=nx*N;
     postype='x';
     nn=[nx*N+2*M 0 2*M*N 0];
     th=[];
  else
     name='RADAR';
     posind=0;
     postype='th';
     nn=[nx*N 0 2*M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2);']];
	     h=[h,['atan2(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,'),x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,'));']];
     end
  end
  h(end)=[];
  h=[h,']'];
  ylabel=repmat({'Range','Bearing'},1,N*M);
  pe=1e-2*eye(nn(3));

elseif strcmp(lower(ex(1:5)),'gps2d');
  name='GPS2d';
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['x(',num2str(2*n-2+1),',:);x(',num2str(2*n-2+2),',:);']];
     end
  end
  h(end)=[];
  h=[h,']'];
  nn=[2*N 0 2*N*M 0];
  ylabel=repmat({'X','Y'},1,M*N);
  xlabel=repmat({'X','Y'},1,N);
  pe=50*eye(nn(3));
  th=[];

elseif strcmp(lower(ex(1:5)),'gps3d');
  name='GPS3d';
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['x(',num2str(3*n-3+1),',:);x(',num2str(3*n-3+2),',:);x(',num2str(3*n-3+3),',:);']];
     end
  end
  h(end)=[];
  h=[h,']'];
  nn=[3*N 0 3*N*M 0];
  ylabel=repmat({'X','Y','Z'},1,M*N);
  xlabel=repmat({'X','Y','Z'},1,N);
  pe=50*eye(nn(3));
  th=[];


elseif strcmp(lower(ex(1:5)),'mag2d');
  if slam
     name='Magnetometer 2D SLAM';
     posind=nx*N;
     postype='x';
     nn=[2 0 3 3];
     th=[];
  else
     name='Magnetometer 2D';
     posind=0;
     postype='th';
     nn=[2 0 3 3];
  end
  mu0=4*pi*1e-7; % Magnetic constant
  th=randn(3,1);
  h=['th(1:3)*1e-7./sqrt(x(1,:).^2+x(2,:).^2).^3.*sqrt(1+3*sin(atan2(x(2,:),x(1,:))).^2)'];
  ylabel={'Bx','By','Bz'};
  xlabel={'X','Y'};

elseif strcmp(lower(ex(1:5)),'mag3d');
    if slam
     name='Magnetometer 3D SLAM';
     posind=nx*N;
     postype='x';
     nn=[3 0 3 3];
     th=[];
  else
     name='Magnetometer 3D';
     posind=0;
     postype='th';
     nn=[3 0 3 3];
  end
  mu0=4*pi*1e-7; % Magnetic constant
  th=randn(3,1);
  h=['1e-7/sqrt(x(1).^2+x(2).^2+x(3).^2)^5*(3*(x(1:3)''*th(1:3))*x(1:3) - (x(1:3)''*x(1:3))*th(1:3))'];
  ylabel={'Bx','By','Bz'};
  xlabel={'X','Y','Z'};

elseif strcmp(lower(ex(1:4)),'quat');
  name='Quaternion constraint';
  if nargin<4
     qind=[1:4];
  else
     qind=varargin{1};
     if ~isvector(qind) | length(qind)~=4
       error('Index argument must be a vector of length 4')
     end
  end
  nx=max([qind(:);4]);
  h=['x(',num2str(qind(1)),',:).^2+x(',num2str(qind(2)),...
      ',:).^2+x(',num2str(qind(3)),',:).^2+x(',num2str(qind(4)),',:).^2-1'];
  nn=[nx 0 1 0];
  th=[];
  x0=zeros(nx,1);
  x0(qind)=rand(4,1);
  x0=x0/norm(x0);
  pe=1e-8;
  ylabel={'0'};
  xlabel=repmat({''},1,nx);
  for j=1:4
    xlabel{qind(j)}=['q',num2str(j)];
  end


elseif strcmp(lower(ex(1:7)),'amundin');
    %   'amundin'  ||pt(n)-ps(m)||+x(N*nx+(m-1)1)+n*x(N*nx+2)
    if slam
     name='Amundin SLAM  with clock offset and drift';
     posind=nx*N+1;
     postype='x';
     nn=[nx*N+2*(M-1) 0 M*2 0];
     th=[];
  else
     name='Amundin with clock offset and drift';
     posind=0;
     postype='th';
     nn=[nx*N+1 0 M*N 2*M];
  end
  h=['['];
  for n=1:N
     for m=1:M
	     h=[h,['sqrt((x(',num2str(2*n-1),',:)-',postype,'(',num2str(posind+2*m-1),posx,')).^2+(x(',num2str(2*n),',:)-',postype,'(',num2str(posind+2*m),posx,')).^2) + x(',num2str(nx*N+1),',:)']];
        if m<M
	   h=[h,';'];
        end
     end
     if n<N
        h=[h,';'];
     end
  end
  h=[h,']'];
  %ylabel={'TOA with offset'};



else
   error(['Unknown sensor model option: ',ex])
end
s=sensormod(h,nn);
if exist('x0')==1
   s.x0=x0;
else
   s.x0=rand(nn(1),1);
end
if exist('th')==1
   s.th=th;
else
   s.th=rand(2*M,1);
end
if exist('pe')==1
   s.pe=pe;
else
   s.pe=1e-4*eye(nn(3));
end
s.name=name;
s.fs=1;  % Default value XXX
if exist('ylabel')==1
   s.ylabel=ylabel;
end
if exist('xlabel')==1
   s.xlabel=xlabel;
end
