function p=confband2(x,y,Py,plottype,col,ax,level,conftype,lw)

% Adds a confidence bound of level (%) to a plot based on covariance
% confband(x,y,Py,plottype,col,ax,level,type)
% conftype 1 gives filled area, 2 gives upper and lower bound

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<9, lw=1; end
if nargin<8, conftype=1; end
if level>1; level=level/100; end  % Per cent or relative?
if level<0 or level>1
   error('level must be in [0,1]')
end
X=ndist(0,1);
c=erfinv(X,level);
yupper=y+c*sqrt(Py);
ylower=y-c*sqrt(Py);
x=x(:);
if conftype==2  % conf lines
    p(1)=feval(plottype,x,yupper,['--',col],'parent',ax,'linewidth',lw);
    p(2)=feval(plottype,x,ylower,['--',col],'parent',ax,'linewidth',lw);
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
    p=patch(xx,yy,cc,'facealpha',0.3,'parent',ax);
    if strcmp(plottype,'semilogy')| strcmp(plottype,'loglog')
        set(ax,'yscale','log')
    end
    if strcmp(plottype,'semilogx')| strcmp(plottype,'loglog')
        set(ax,'xscale','log')
    end
end
