function setpaths(~)
% SETPATHS  Set paths for TSRT14
%
% Set the paths needed in the course TSRT14 to the MATLAB path.  The
% added paths contain the SigSys toolbox and some mat-files used in the
% course.
%
% If given with an argument, remove the paths again.

   cwd = fileparts(mfilename('fullpath'));
 % cwd = '/courses/tsrt14';
 
  %cwd = '/Desktop/TSRT14/sigsys28-Oct-2019open';
   %cwd = ''; 
  if nargin == 0
    add2path(cwd);
    add2path(fullfile(cwd, 'classes'));
    add2path(fullfile(cwd, 'mfiles'));
    add2path(fullfile(cwd, 'data'));
  else
    removepath(fullfile(cwd, 'data'));
    removepath(fullfile(cwd, 'mfiles'));
    removepath(fullfile(cwd, 'classes'));
    removepath(cwd);
  end
end

function add2path(p, doRecurse)
  if nargin<2 || ~doRecurse
    fprintf('Add path: %s\n', p);
    addpath(p);
  else
    fprintf('Recursively add path: %s\n', p);
    addpath(genpath(p));
  end
end

function removepath(p, doRecurse)
  if nargin<2 || ~doRecurse
    fprintf('Remove path: %s\n', p);
    rmpath(p);
  else
    fprintf('Recursively remove path: %s\n', p);
    rmpath(genpath(p));
  end
end
