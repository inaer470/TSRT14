function ok=isnumericscalar(x)
%ISNUMERICSCALAR tests if the argument is numeric and scalar
%
% ok=isnumscalar(x)

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

ok = isnumeric(x) && isscalar(x);
