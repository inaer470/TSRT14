function y=xbow2sig(file)
%xbow2sig converts a log file from the xbow app xsensor to a sig object

if ~exist(file,'file')
    error([file,' not found'])
end
fid=fopen(file,'r');
a=fread(fid);
fclose(fid);
astr=char(a');
if ~strcmp(astr(1:6),'iPhone')
    disp('Not a log file from Crossbow')
end

ind=find(a==10); % carriage return
desc=astr(1:ind(1)-1);
for k=3:length(ind)
    ykstr=astr(ind(k-1)+1:ind(k)-1);
    indnan=findstr(ykstr,'(null)');
    for i=length(indnan):-1:1
        ykstr(indnan(i):indnan(i)+5)='NaN   ';
    end
    x(k-2,:)=str2num(ykstr);
end
t=x(:,1);
order=[5 6 7 11 12 13 2 3 4 8 9 10 14 15 16];
y=x(:,order);
ylab={'a_x [g]','a_y [g]','a_z [g]',...
      '\omega_x [rad/s]',['\omega_y [rad/s]'],'\omega_z [rad/s]',...
      'm_x [muT]','m_y [muT]','m_z [muT]',...
      'lat [deg]','long [deg]','alt [deg]',...
      '\theta [rad]','\phi [rad]','\psi [rad]'};
y=sig(y,t);
y.desc=desc;
y.ylabel=ylab;
fs=1/(t(2)-t(1));
y.tlabel=['Time [s] (f_s=',fs,')'];

