function p=confband(x,yMC,plottype,col,ax,level,type,lw)

% Adds a confidence bound of level (%) to a plot based on MC data
% confband(x,yMC,plottype,col,ax,level,type)
% type 1 gives filled area, 2 gives upper and lower bound

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<8, lw=1; end
if nargin<7, type=1; end
if level<1; level=100*level; end  % Per cent or relative?
MC=size(yMC,2);
yMC=sort(yMC,2);
yupper=yMC(:,round(MC*level/100));
ymedian=yMC(:,round(MC/2));
ylower=yMC(:,ceil(MC*(1-level/100)));
x=x(:);
if type==2  % conf lines
%    p(1)=plot(x,yupper,['--',col],'parent',ax);
%    p(2)=plot(x,ylower,['--',col],'parent',ax);
%    p(3)=plot(x,ymedian,['-.',col],'parent',ax);
    p(1)=feval(plottype,x,yupper,['--',col],'parent',ax,'linewidth',lw);
    p(2)=feval(plottype,x,ylower,['--',col],'parent',ax,'linewidth',lw);
    p(3)=feval(plottype,x,ymedian,['-.',col],'parent',ax,'linewidth',lw);
else  % conf band
    xx=[x(1:end-1) x(2:end) x(2:end) x(1:end-1)]';
    yy=[yupper(1:end-1) yupper(2:end) ylower(2:end) ylower(1:end-1)]';
    %size(xx),size(yy)
     if strcmp(col,'b'), cc=[0.4 0.4 1]; end
     if strcmp(col,'g'), cc=[0.4 1 0.4]; end
     if strcmp(col,'r'), cc=[1 0.4 0.4]; end
     if strcmp(col,'m'), cc=[1 0.4 0.4]; end
     if strcmp(col,'y'), cc=[1 0.4 1]; end
     if strcmp(col,'c'), cc=[0.4 1 1]; end
     if strcmp(col,'k'), cc=[0.4 0.4 0.4]; end
     p=patch(xx,yy,cc,'facealpha',0.3,'parent',ax,'edgecolor','none');
    if strcmp(plottype,'semilogy')| strcmp(plottype,'loglog')
        set(ax,'yscale','log')
    end
    if strcmp(plottype,'semilogx')| strcmp(plottype,'loglog')
        set(ax,'xscale','log')
    end
end
