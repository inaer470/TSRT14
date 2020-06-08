function w=getwindow(N,type,n)

%GETWINDOW computes data windows
%   w=getwindow(N,type,n)   window design syntax
%   N length of window
%   type is 'box', 'hanning', 'hamming','kaiser',
%           'blackman', 'bartlett' or 'spline'
%   spline uses a uniform window convolved with itself n times.
%
%   Examples:
%     w1=getwindow(40,'hanning');
%     w2=getwindow(40,'hamming');
%     subplot(2,1,1), plot([w1 w2])
%     W1=ft(sig(w1)); W2=ft(sig(w2));
%     subplot(2,1,2),
%     semilogy(W1,W2,'type',3,'Ylim',[1e-3 40])
%   See also: sig.window

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2
  type='hanning';
end

switch lower(type)
case 'box'
  w = ones(N,1);
case 'hanning'
  w = 0.5*(1 - cos(2*pi*(1:N)'/(N+1)));
case 'hamming'
  a=0.54;
  w = a - (1-a)*cos(2*pi*(1:N)'/(N+1));
case 'kaiser'
    if nargin<3;
        alpha=3;
    else
        alpha=n;
    end
    w=besseli(0,alpha*sqrt(1-(2*(1:N)/(N+1)-1).^2))/besseli(0,alpha);
    w=w';
case 'blackman'
  w = 0.42-0.5*cos(2*pi*(1:N)'/(N+1))+0.08*cos(4*pi*(1:N)'/(N+1));
case 'bartlett'
  if N/2==round(N/2)
    w=[2*(1:N/2)'/(N+1); 2-2*(N/2+1:N)'/(N+1)];
  else
    w=[2*(1:(N+1)/2)'/(N+1); 2-2*((N+1)/2+1:N)'/(N+1)];
  end
case 'spline'
    if nargin<3;
        n=3;
    end
  m=N/2^n;
  if m~=round(m)
       error('N must be divisible with 2^n')
  else
       wtmp=ones(m,1);
       w=wtmp;
       for i=2:n
           w=[conv(w,wtmp); 0];
       end
       w=w/max(w);
  end
end
