function s=mat2strformat(n,format);

%MAT2STR is the true inverse of str2num which works for matrices
%in contrast to num2str.
%   s=mat2str(n,format)
% Optional format is the same as in num2str {'%11.4g'}

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<2; format='%11.4g'; end
[nr,nc]=size(n);

if nr>1 | nc>1
  s='[';
  for i=1:nr
    for j=1:nc
      s=[s,num2str(n(i,j),format),','];
    end
    s(length(s))=';';
    s=[s,''];
  end
  s(length(s))=']';
else
  s=num2str(n,format);
end
