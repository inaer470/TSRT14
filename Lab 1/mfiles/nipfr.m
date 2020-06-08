function dx=nipfr(t,x,u,th)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

X=x(1,:);
T=x(2,:);

k1=exp(th(7) + th(8)./T);
k2=exp(th(9) + th(10)./T);
Keq=th(11)+th(12)./T;
Rn=X.*sqrt(1-th(2).*(1-X)) - th(3)*(1-X)/Keq;
Rd=(k1+k2*(1-X)).^2;
R=Rn/Rd;

dx1=th(1)*R;
dx2=th(4)*(T-th(5)) + th(6)*R;
dx=[dx1;dx2];
