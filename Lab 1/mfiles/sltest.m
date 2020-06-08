function sltest

%close all
load sltest
N=size(y);
z=y.y;
z=z-ones(N,1)*mean(z,1);
T=1/32;
g=[0; 0; -9.81];

ang=[0;0;0];
P=ones(3);
Q=1e-3*eye(3);
R=1e-3*eye(3);

Ang=zeros(N,3);
for t=60:N
    % time update
    w=y.y(t,4:6)';%*pi/180;
    w(1)=-w(1);
    F=[1 -sin(ang(1))*tan(ang(2)) cos(ang(1))*tan(ang(2));...
       0 cos(ang(1)) sin(ang(1));...
       0 -sin(ang(1))*sec(ang(2)) cos(ang(1))*sec(ang(2))];
    ang=ang+T*F*w;
    P=P+F*Q*F';  % gyro noise 
    % measurement update
    a=9.81*y.y(t,1:3)';
    a(3)=-a(3);
    %    g, Qrot(ang,a),a,pausef
    epsi=g-Qrot(ang,a);
    H=Jrot(ang,a);
    K=P*H'*inv(H*P*H'+R);
    ang=ang+K*epsi;
    P=P-K*H*P;
    Ang(t+1,:)=ang';
    ang0=y.y(t,13:15)*180/pi;
    anghat=ang'*180/pi;
    %    ang0, anghat, pause
    
end
figure(1)
plotfix
plot(Ang)

return
L=20:50;
ca=zeros(N, L(end));
cg=zeros(N, L(end));
for t=L(end)+1:N-L(end);
    for l=L
        for i=1:3;
           ca(t,l)=ca(t,l)+z(t+1:t+l,i)'*z(t-l+1:t,i) / ...
                  (z(t-l+1:t+l,i)'*z(t-l+1:t+l,i));
        end
        for i=4:6;
           cg(t,l)=cg(t,l)+z(t+1:t+l,i)'*z(t-l+1:t,i) / ...
                  (z(t-l+1:t+l,i)'*z(t-l+1:t+l,i));
        end
    end
end
[camax,indamax]=max(ca');
[cgmax,indgmax]=max(cg');
figure(2)
plotfix
plot((1:N)/32,[indamax' indgmax']/32)
legend('acc','gyro')
xlabel('Time')
ylabel('Step cycle time')

end

function h=Qrot(ang,a)
    Qz=[1 0 0; 0 cos(ang(1)) sin(ang(1));  0 -sin(ang(1)) cos(ang(1))];
    Qy=[cos(ang(2)) 0 -sin(ang(2));  0 1 0; sin(ang(2)) 0  cos(ang(2))];
    Qx=[cos(ang(3)) sin(ang(3)) 0; -sin(ang(3)) cos(ang(3)) 0; 0 0 1];
    h=Qz*Qy*Qx*a;
end

function H=Jrot(ang,a)
    s=1e-6;
    H0=Qrot(ang,a);
    H1=Qrot(ang+s*[1;0;0],a);
    H2=Qrot(ang+s*[0;1;0],a);
    H3=Qrot(ang+s*[0;0;1],a);
    H=[H1-H0 H2-H0 H3-H0]/s;
end
    
    
    
    