function y = sniper1(t,x,u,th)
% y = sniper1(t,x,u,th)
% Calculates the acoustic event times for the
% muzzle blast y(1:M) and shock wave y(M+1:2M)
% for M microphones at positions th (M-by-2).
% The shooter is at position s and fires at time
% t0 in the direction a. Bullet speed is v, speed
% of sound is c.
%
% The parameter vector is
% x = [s(1) s(2) t0 c v a]
%
% y(M+k) = NaN if microphone k is not reached by
% the shock wave.
%
% David Lindgren 090319,090320

  % Bullet deceleration constant
  r = 0.63; % Must be >0 !

 s(1) = x(1);     % Shooter pos east
 s(2) = x(2);     % Shooter pos north
   t0 = x(3);     % Firing time
    c = x(4);     % Speed of sound
   vb = x(5);     % Bullet speed
    a = x(6);     % Shooting angle

  M = length(th)/2;
  th=reshape(th,2,M)';

  % Microphone-shooter vectors
  z = th-repmat(s,M,1);

  % Muzzle blast
  yMB = t0 + sqrt(sum(z.^2, 2))/c;

  % Shock wave
  gamma = abs(atan2(z(:,2),z(:,1)) - a);
  for i = 1:M,
    d = swGenPt(gamma(i), norm(z(i,:)), vb, r, c);
    if ~isnan(d),
      ySW(i,1) =  t0 +...
          log(vb/(vb-r*d))/r  +...
          norm(s+[cos(a) sin(a)]*d-th(i,:))/c;
    else
      ySW(i,1) = NaN;
    end
  end
  y = 340*[yMB;ySW];


function dk = swGenPt(a, p, v, k, c),
% Finds the distance from shooter to the
% point where the shockwave that reaches
% p has its origin.
%
% a Angle between LOS shooter/sensor and
%   shooting direction.
% p Microphone/shooter distance.
% v Initial bullet speed.
% k Deceleration parameter (=r).
% c Speed of sound.

  if v<=c, % No SW
    dk = NaN;
    return
  end

  beta = asin(c/v);
  if a>(90-beta), % SW not within cone
    dk = NaN;
    return
  end

  g = 1/p/sin(a);
  e = -cos(a)/sin(a);
  f =  [g^2*k^2
        -2*g^2*v*k+2*g*e*k^2
        g^2*v^2+e^2*k^2-4*g*e*v*k-g^2*c^2
        2*g*e*v^2-2*e^2*v*k-2*c^2*g*e
        e^2*v^2-c^2-c^2*e^2];
  R= roots(f);
  R = R(~imag(R));

  % Check if reasonable solution, else not at SW

  if isempty(R), % Noo real roots
    dk = NaN;
    return;
  end
  dk = min(R); % Least real root reasonable

  % Check v>c at dk
  if v-dk*k<=c || dk<=0,
    dk = NaN;
  end
