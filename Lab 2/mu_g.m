function [x,P] = mu_g(x, P, yacc, Ra, g0)
% Implements the accelerometer measurement updates function

[h_derv1, h_derv2, h_derv3, h_derv4] = dQqdq(x);
h_derv = [h_derv1*g0, h_derv2*g0, h_derv3*g0, h_derv4*g0]; % 13.20 d
Q = Qq(x); % 13.16

S = Ra + h_derv*P*h_derv'; % 8.3a
K = P*h_derv'*inv(S); % 8.3b
epsil = yacc - Q*g0; % 8.3c   13.19a

x = x + K*epsil; % 8.3d
P = P-P*h_derv'*inv(S)*h_derv*P; % 8.3e

[x, P ] = mu_normalizeQ(x, P); % quaternion keeps unit length

end

