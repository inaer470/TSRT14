function sout=inherit(sout,s,nameext)
%INHERIT transfers public fields from one object to another
%   sout=inherit(sout,s,nameext)
%   nameext is extended to the object name to trace its background

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

if nargin<3, nameext=[]; end

sout.MC=s.MC;
sout.name=s.name;
% XXX sout.desc=char(s.desc,nameext);
try
  if ~isempty(s.ulabel)
    sout.ulabel=s.ulabel;
  end
end
%s.ylabel,sout.ylabel
try
  if ~isempty(s.ylabel)
    sout.ylabel=s.ylabel;
  end
end
try
  if ~isempty(s.xlabel)
    sout.xlabel=s.xlabel;
  end
end
try
  if ~isempty(s.tlabel)
    sout.xlabel=s.tlabel;
  end
end
try
  if ~isempty(s.marker)
    sout.marker=s.marker;
  end
end
try
  if ~isempty(s.markerlabel)
    sout.xlabel=s.markerlabel;
  end
end
