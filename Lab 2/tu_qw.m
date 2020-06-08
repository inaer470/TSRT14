function [x, P] = tu_qw(x, P, omega, T, Rw)
% Implements the time update function
% EKF TT1, p 197

Q = T^2/4*Sq(x)*Rw*Sq(x)';

if nargin < 5
    
    P = Q + P;
    
else
    I = eye(length(x));
    f = (I + 1/2*Somega(omega)*T); % eq. 2k
    x = f * x;
    
    P = Q+f*P*f';
end

[x, P] = mu_normalizeQ(x, P); % quaternion keeps unit length

end

