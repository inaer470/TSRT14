function [p,fs]=randpoly(n,varargin);

%RANDPOLY generates a random polynomial of specified order
%   [p,fs]=randpoly(n,varargin);
%
%   The polynomial has all roots inside the unit circle.
%   The phase of the roots is uniformly distributed (avoiding
%   the negative real axis), while the radius
%   is in the interval [0,0.9] and [0,0.99], respectively, with most of the
%   distribution closer to the unit circle than to 0.
%   This yields high variation in the spectra with limited peaks.
%
%   With fs=0, the roots are distributed in the left half plane.
%
%
%   Property  Value   Description
%   frac      {0.7}   Fraction of real poles and zeros
%   fs        {2}     Sampling frequency
%
%   The function is called by rand in the ss and tf objects

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

opt=struct('frac',0.7,'fs',NaN);
opt=optset(opt,varargin);

fs=opt.fs;
% Number of complex poles even and random between 0 and n
nc=min([2*floor(n/2) max([0 2*round(opt.frac*n+n/3*randn(1))])]);
nr=n-nc; % real poles
if fs>0 %Discrete time, poles and zeros inside unit circle
    pr=rand(1,nr);
    pc=0.9*(1-rand(1,nc/2).^2).*exp(i*pi*rand(1,nc/2));
    pc=[pc conj(pc)];
    poles=[pr pc];
    p=real(poly(poles));
else
    pr=-log(rand(1,nr)+0.01);
    pc=-log(rand(1,nc/2)+0.01);
    phasec=0.9*rand(1,nc/2)*pi/2;
    phasec=[pi+phasec pi-phasec];
    poles=[-pr [pc pc].*exp(i*phasec)];
    p=real(poly(poles));
end
