function x = ball(t,x,u,th)
%BALL gives the dynamics of a ball subject to gravitation and the Magnus effect
%  nn=[9 0 ny 4];
%  th=[th1 th2 th3 T],
%  th1 screw coefficient,
%  th2 air drag
%  th3 damping of angular speed
%  x=[p,v,omega]
%  u=[]
%  y=[x(2);x(3)] as an example
%
% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


g=[0;0;9.81];
T=th(4); th1=th(1); th2=th(2); th3=th(3);
v=x(4:6);
omega=x(7:9);
S = [0 -omega(3) omega(2);
     -omega(3) 0 omega(1);
     -omega(2) omega(1) 0];
x = x + T*[v;-th(1)*(v'*v)*v+th2*S*v-g;-th3*omega];
end
