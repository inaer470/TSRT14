function [x,P] = mu_m(x, P, ymag, Rm, m0)
% Implements magnetometer measurement update
[h_derv1, h_derv2, h_derv3, h_derv4] = dQqdq(x);
h_derv = [h_derv1*m0, h_derv2*m0, h_derv3*m0, h_derv4*m0];
Q = Qq(x);

S = Rm + h_derv*P*h_derv'; % 8.3a
K = P*h_derv'*inv(S); % 8.3b
epsil = ymag - Q*m0; % 8.3c   13.19a

x = x + K*epsil; % 8.3d
P = P-P*h_derv'*inv(S)*h_derv*P; % 8.3e

[x, P ] = mu_normalizeQ(x, P); % quaternion keeps unit length

end

