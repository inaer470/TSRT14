function fp=numgrad(f,ind,varargin)
%NUMGRAD returns a numeric gradient approximation
%   fp=numgrad(f,ind,varargin)
%   Central difference approximation
%         fp(x)=(f(x+h)-f(x+h))/2h
%   of f as a function of varargin{ind}
%
%   Normally, ind=1, and the gradient with respect to the first element
%   is returned. That is, if f=f(x), then fp=df(x)/dx
%   For a function as f(t,x,u,th) (standard form in the NL object)
%   fp=numgrad(f,2,t,x,u,th) returns fp=df(x)/dx
%
%   Example:
%   x=1; f='[x^2+x+x^3+x^4]',numgrad(f,1,x)
%   x=[1 1]'; y=2; f='[x''*x+y*x(2)^10]',numgrad(f,1,x,y)
%   x=[1 1]'; y=2; f='[x''*x+y^2*x(2)^10]',numgrad(f,2,x,y)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


arg=varargin;
if isstr(f) & exist(f)~=2
    f=inline(f);
end
try
   fpx=feval(f,arg{:});
catch
   error('NUMGRAD: f(varargin{:}) cannot be evaluated')
end
h=max([sqrt(eps)*norm(arg{ind}) sqrt(eps)]);
n=length(arg{ind});
I=eye(n);
arg1=arg;
arg2=arg;
for i=1:n
   arg1{ind}=arg{ind}+h*I(:,i);
   arg2{ind}=arg{ind}-h*I(:,i);
   fp(:,i)=(feval(f,arg1{:})-feval(f,arg2{:}))/2/h;
end
