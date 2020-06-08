function out=mchar(p,var,format,timesstr)
%MCHAR generates polynomial strings

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if nargin<3, format='%11.4g'; end
if nargin<4, timesstr='*'; end

if all(p==0)
  out = '0';
else
  d = length(p)-1;
  out = [];
  for a=p
    if a~=0
      if ~isempty(out)
        if a>0
          out = [out '+'];
        else
          out = [out '-'];
          a = -a;
        end
      end
      if a~=1 | d==0
        out = [out num2str(a,format)];
        if d>0
          out = [out timesstr];
        end
      end
      if d>=2
        out = [out,var,'^' int2str(d)];
      elseif d==1
        out = [out,var];
      end
    end
    d = d-1;
  end
end
if ~(isempty(findstr(out,'+')) & isempty(findstr(out,'-')) )
%    out=['(',out,')'];
end
