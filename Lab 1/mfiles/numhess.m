function fpp=numhess(f,ind,varargin)
%NUMHESS returns a numeric Hessian approximation
%   fp=numhess(f,ind,varargin)
%   Central difference approximation
%   fpp(i,j)=(f(x+h*ei+h*ej)-f(x+h*ei-h*ej))-(f(x-h*ei+h*ej)-f(x-h*ei-h*ej))/4h^2
%   of f as a function of varargin{ind}
%
%   Normally, ind=1, and the hessient with respect to the first element
%   is returned. That is, if f=f(x), then fpp=d2f(x)/dx2
%   For a function as f(t,x,u,th) (standard form in the NL object)
%   fp=numhess(f,2,t,x,u,th) returns fp=df(x)/dx
%
%   Example:
%   x=1; f='[x^2+x+x^3+x^4]',numhess(f,1,x)
%   x=[1 1]'; y=2; f='[x''*x+y*x(2)^10]',numhess(f,1,x,y)
%   x=[1 1]'; y=2; f='[x''*x+y^2*x(2)^10]',numhess(f,2,x,y)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

arg=varargin;

if isstr(f), f=inline(f); end
try
   fpx=feval(f,arg{:});
catch
   error('NUMHESS: f(varargin{:}) cannot be evaluated')
end

h=max([sqrt(sqrt(eps))*norm(arg{ind}) sqrt(sqrt(eps))]);
n=length(arg{ind});
I=eye(n);
arg1=arg;
arg2=arg;
arg3=arg;
arg4=arg;
for i=1:n
   for j=1:i
      arg1{ind}=arg{ind}+h*I(:,i)+h*I(:,j);
      arg2{ind}=arg{ind}+h*I(:,i)-h*I(:,j);
      arg3{ind}=arg{ind}-h*I(:,i)+h*I(:,j);
      arg4{ind}=arg{ind}-h*I(:,i)-h*I(:,j);
      f1=(feval(f,arg1{:})-feval(f,arg2{:}))/2/h;
      f2=(feval(f,arg3{:})-feval(f,arg4{:}))/2/h;
      fpp(:,i,j)=(f1-f2)/2/h;
      fpp(:,j,i)=fpp(:,i,j);
   end
end
