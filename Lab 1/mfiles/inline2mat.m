function c=inline2mat(g,ginit,gend)
%INLINE2MAT creates a character matrix from an inline object
%   c=inline2mat(g,ginit,gend)
%   Symbolically, c=[ginit,g,gend]
%   but the math rows in g are expanded to char matrix

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

g=formula(g);
ind1=find(g=='[');
ind2=find(g==']');
if strfind(g,']*x'), gx='*x'; else, gx=[]; end
g=g(ind1(1)+1:ind2(end)-1);
ind3=find(g==';');
gstr=g(1:ind3(1)-1);
if length(ind3)>0
  for k=2:length(ind3)
	gstr=char(gstr,g(ind3(k-1)+1:ind3(k)-1));
  end
  k=length(ind3);
  if length(ind3)>0
     gstr=char(gstr,g(ind3(k)+1:end));
  end
  [n,m]=size(gstr);
  if n>1
    g1=['/ ';repmat('| ',n-2,1);'\ '];
    g2=[' \';repmat(' |',n-2,1);' /'];
    c=[g1 gstr g2 gx];
  end
  if nargin>1  & ~isempty(ginit);
    n1=size(ginit,1);
    cinit=char(repmat(' ',ceil((n-n1)/2),1),ginit,repmat(' ',floor((n-n1)/2),1));
    if n-n1==1
       cinit(end,:)=[];
    end
    c=[cinit,c];
  end
  if nargin>2 & ~isempty(gend);
     n1=size(gend,1);
     cend=char(repmat(' ',ceil((n-n1)/2),1),[gx,gend],repmat(' ',floor((n-n1)/2),1));
     if n-n1==1
        cend(end,:)=[];
     end
     c=[c cend];
  end
end
