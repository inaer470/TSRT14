function [a,b,c]=thnn2abc(th,nn,P)

%THNN2ABC converts parameter and index vectors to a, b and c polynomials
%   [a,b,c]=thnn2abc(th,nn,P)
%   Internal function for arx and armax objects

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<3; P=[]; end
nn=nn(:)';
th=th(:);
if nargout==2 % arx form
    nn=[nn(1:2) 0 nn(3:end)];
end
if length(nn)<6  % SISO
    nn=[nn 1 1];
end

if P==0; P=[]; end
if ~size(nn,2)==4
     error('Model structure nn must have four elements')
end
if ~size(th,1)==sum(nn(1:3))
     error('Number of rows in th must be the sum of nn(1:3)')
end
na=nn(1); nb=nn(2); nc=nn(3); nk=nn(4); nu=nn(5); ny=nn(6);

a=[1 th(1:na).' ];
if nb==0
    b=[];
else
    for k=1:ny
       b(k,nk+1:nk+nb,:)=reshape(th(na+1+(k-1)*nu*nb:na+nb*nu*k),nb,nu);
    end
end
if nc<1
     c=1;
else
    for k=1:ny
       c(k,:,:)=reshape(th(na+1+nu*nb*ny+(k-1)*ny*nc:na+nu*nb*ny+nc*ny*k),nc,ny);
    end
    cij=th(na+nc*ny*ny+1:end); %'
    c=[ones(ny,1) reshape(cij,ny,ny).'];    %'
end
