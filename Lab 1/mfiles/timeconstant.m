function [Ts,Tf]=timeconstant(a,fs)

%TIMECONSTANT approximates the time constants based on the dominating poles
%   [Ts,Tf]=timeconstant(a)
%   Ts is the slowest time constant
%   Tf is the slowest time constant
%   a is either the denominator polynomial of a TF or A matrix in a SS model

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2
   error('TIMECONSTANT: usage [Ts,Tf]=timeconstant(a,fs)')
end
if size(a,1)==size(a,2)
    r=eig(a);
elseif size(a,1)==1
    r=roots(a);
else
    error('TIMECONSTANT: a is either a row vector or square matrix')
end
if isnan(fs); % Continuous time system
   ind=find(real(r)<0);  % Remove integrators and unstable poles
   if ~isempty(ind)
      r=r(ind);
      Ts=30/min(abs(r));% Dominating pole
      Tf=30/max(abs(r));% Dominating pole
   else
      Ts=10; Tf=10;
   end
else
   r=abs(r);
   ind=find(r<1);   % Remove integrators and unstable poles
   if ~isempty(ind)
       Ts=log(100)/min(1-r(ind)); % In samples
       Tf=log(100)/max(1-r(ind));
   else
       Ts=30; Tf=30;
   end
   Ts=Ts/fs;
   Tf=Tf/fs;
end
