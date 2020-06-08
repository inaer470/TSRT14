function [opt,errmsg]=optset(opt,varargin);

%OPTSET overwrites default values in opt by user supplied values
%   [opt,errmsg]=optset(opt,Property1,Value1,...);
%
%   errmsg contains error messages

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


if isfield(opt,'fontsize')
    global SIGNAL_FONTSIZE
    opt.fontsize=SIGNAL_FONTSIZE;
    if isempty(opt.fontsize)
        opt.fontsize=12;
    end
end
if isfield(opt,'linewidth')
    global SIGNAL_LINEWIDTH
    opt.linewidth=SIGNAL_LINEWIDTH;
    if isempty(opt.linewidth)
        opt.linewidth=1;
    end
end
if isfield(opt,'markersize')
    global SIGNAL_MARKERSIZE
    opt.markersize=SIGNAL_MARKERSIZE;
    if isempty(opt.markersize)
        opt.markersize=5;
    end
end
varargin=varargin{1};
n=length(varargin)/2;
if n~=round(n)
    disp(['Entered parameter value pairs: ',varargin{:}])
    error('Number of parameter and values for options must be even')
end
for k=1:n
    if ~isstr(varargin{2*k-1})
        error('Optional arguments must be ordered in property (string) and value pairs')
    end
    if ~isfield(opt,varargin{2*k-1})
        errmsg=['Unknown option: ',varargin{2*k-1}];
        if nargout<2; error(errmsg), end
    else
        opt=setfield(opt,varargin{2*k-1},varargin{2*k});
    end
end
