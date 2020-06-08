function y=rsssls(t,x,u,th)
%RSSSLS implements the separable RSS measurement model
%   th positions of the sensors
%   x  positions of the targets
%   u  =y, used for estimation of nuisance parameters

M=length(th)/2;
K=length(x)/2;

for k=1:K;
    xk=x(k,:)';
    for m=1:M;
       pm=th(m*2-1:2*m);
       c(m,1)=0.5*log((xk(1)-pm(1))^2+(xk(2)-pm(2))^2);
    end
    phi=[ones(M,1) c];
    R=phi'*phi;
    f=phi'*u;
    Rinv=inv(R);
    thhat=Rinv*f;
% unfinished
end
