function plotfix(h,fs,lw)

%PLOTFIX
%   plotfix(h,fs,lw)
%   Change the default font size and line width in figure.
%   h   gcf      handle, 0 makes the definition
%                global to all new figures.
%   fs 16        fontsize
%   lw 2         linewidth

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if nargin<3, lw=2; end
if nargin<2, fs=16; end
if nargin<1, h=gcf; end
if length(lw)==1; lw=[lw lw]; end
if length(fs)==1; fs=[fs fs]; end

set(h,...
 'DefaultAxesFontSize',fs(1),...
 'DefaultLineLineWidth',lw(1),...
 'DefaultAxesLineWidth',lw(2),...
 'DefaultAxesGridLineStyle',...
      '--',...
 'DefaultTextFontSize',fs(2));
