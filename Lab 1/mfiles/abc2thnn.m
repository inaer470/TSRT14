function [th,nn,P]=abc2thnn(a,b,c)
%ABC2THNN converts a,b and c polynomials to parameter and index vectors
%   [th,nn,P]=abc2thnn(a,b,c)
%   Internal function used in arx and armax objects
%

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


na=size(a,2)-1;
nb=size(b,2);
if nb>0
   ny=size(b,1);
   nu=size(b,3);
else
   ny=1;
   nu=0;
end
%nc=size(c,2);
MC=size(b,4)-1;

th=[];
P=[];
b=b/a(1);
a=a/a(1);
if na>0
    th=[a(1,2:end)].';
end

if nb>0
    for i=1:nu
        for j=1:ny
            ind=find(b(j,:,i)~=0);
            NK(j,i)=ind(1)-1;
        end
    end
    nk=min(min(NK));
    nb=nb-nk;
    for k=1:ny
       bb=b(k,nk+1:end,:);
       th(na+1+(k-1)*nu*nb:na+nb*nu*k)=bb(:);
    end
else
    nk=0;
end

if isempty(c)
    nc=[];
else
    for k=1:ny
       cc=c(k,2:end,:);
       th(na+1+nu*nb*ny+(k-1)*ny*nc:na+nu*nb*ny+nc*ny*k)=cc(:);
    end
end
nn=[na nb nc nk nu ny];
th=th(:);
P=zeros(length(th));
