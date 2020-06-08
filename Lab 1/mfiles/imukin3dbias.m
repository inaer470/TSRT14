function x = imukin3dbias(t,x,u,th)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

  % Kinematic 3D model for IMU input with six bias states
  T=th(1);
  a=u(1:3)+x(11:13);
  omega=u(4:6)+x(14:16);
  F = [eye(3)     T*eye(3)   zeros(3,4);
        zeros(3)   eye(3)     zeros(3,4);
        zeros(4,3) zeros(4,3) Somega(omega,T)];
  Gu=[.5*T^2*eye(3); T*eye(3); zeros(4,3)];
  x(1:10) = F*x(1:10) + Gu*(Rq(x(7:10)).'*a - ginert(x(1:3)')');
end

function S=Somega(omega, T)
% The result of this function is the system matrix for the quaternion
% dynamics
  c1 = cos(norm(omega)*T/2);
  c2 = 1/2*sinc(norm(omega)*T/2);
  crossom = [0 omega(3) -omega(2); -omega(3) 0 omega(1); omega(2) -omega(1) 0];
  S = (c1*eye(4) + c2*[0 -omega.'; omega crossom]);
end

function R=Rq(q)
   q0=q(1);   q1=q(2);   q2=q(3);   q3=q(4);
   R = [(q0^2+q1^2-q2^2-q3^2) 2*(q1*q2+q0*q3)       2*(q1*q3-q0*q2);
         2*(q1*q2-q0*q3)      (q0^2-q1^2+q2^2-q3^2) 2*(q2*q3+q0*q1);
         2*(q1*q3+q0*q2)      2*(q2*q3-q0*q1)       (q0^2-q1^2-q2^2+q3^2)];
end

function gx = gravx (xi)
  if all(xi==0)
      gx=zeros(3,1);  % special case for NL constructor
      return
  end
  mu = 3.986013e14;   % universal gravity constant
  J2 = 0.00108263;    % J2 gravity constant
  Re = 6378137.0;     % Earth's equator radius
  ri = 1/norm(xi);
  xir = xi*ri;
  xx = xir(:)*xir(:)';
  xx5 = 0.5*xx;
  xx5(:,3) = 5*xx5(:,3);
  gxJ2 = 5*(xx(3,3)*(7*xx - eye(3)) - xx5 - xx5') + diag([1 1 3]);
  gx = mu*ri^3*(eye(3) - 3*xx + 1.5*J2*(Re*ri)^2*gxJ2);
end

function g=ginert(inert)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants and parameters
Re = 6378137.0;            % Earth's equator radius
mu = 3.986013e14;          % universal gravity constant
J2 = 0.00108263;           % J2 gravity constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = sqrt(sum(inert.^2,2));

x=inert(:,1);
y=inert(:,2);
z=inert(:,3);
g(:,1)=(1/2)*mu*x.*(-12*J2*Re^2*z.^2 + 3*J2*Re^2*x.^2 + 3*J2*Re^2*y.^2 + ...
   2*x.^4 + 4*x.^2.*y.^2 + 4*x.^2.*z.^2 + 2*y.^4 + 4*y.^2.*z.^2 + 2*z.^4)./r.^7; g(:,2)=(1/2)*mu*y.*(-12*J2*Re^2*z.^2+3*J2*Re^2*x.^2+3*J2*Re^2*y.^2+2*x.^4+4*x.^2.*y.^2+4*x.^2.*z.^2+2*y.^4+4*y.^2.*z.^2+2*z.^4)./r.^7;
g(:,3)=(1/2)*mu*z.*(-6*J2*Re^2*z.^2 + 9*J2*Re^2*x.^2 + 9*J2*Re^2*y.^2 + ...
   2*x.^4 + 4*x.^2.*y.^2 + 4*x.^2.*z.^2 + 2*y.^4 + 4*y.^2.*z.^2 + 2*z.^4)./r.^7;
end
