function texcode=texmatrix(V,varargin)

%Generate a LaTex table from a matlab matrix
%   texcode=texmatrix(V,Property1,Value1,...)
%
%   V is the matrix to be coded.
%   The output is the produced texcode that can be pasted into
%   any latex document. Alternatively, the code is put into a file,
%   which is inputted by reference into the document.
%   Filters (freeware or shareware) for other word processors are available:
%   - texpoint: freeware for Powerpoint
%   - LaImport: FrameMaker
%   - tex2word: Word
%   - latex2rtf: RTF documents
%   - tex4ht: HTML or XML hypertext documents
%
%   Property   Value/{Default}  Description
%   --------------------------------------------------------------
%   filename   {''}             Name of the .tex file (none for '')
%   decimals   {1}              Number of decimals
%   env        {'eqnarray*'}    Tex environment, '' means no env
%
%   Examples:
%
%   See also: textable, ss/tex, tf/tex

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

opt=struct('filename','','decimals',1,'env','eqnarray*');
opt=optset(opt,varargin);

[nrow,ncol]=size(V);

if ~isempty(opt.env)
    texcode=sprintf(['\\begin{',opt.env,'}']);
else
    texcode='';
end
if nrow>1 | ncol>1
    texcode=strvcat(texcode,sprintf('\\left['));
end
tabarg='r';
for i=2:ncol;
  tabarg=[tabarg,'r'];
end;
texcode=strvcat(texcode,sprintf(['\\begin{array}{',tabarg,'}']));
for j=1:nrow;
  rowj=['%.',num2str(opt.decimals),'f'];
  for i=2:ncol;
    rowj=[rowj,' & %.',num2str(opt.decimals),'f'];
  end;
  if j<nrow; rowj=[rowj,' \\\\']; end
  texcode=strvcat(texcode,sprintf(rowj,V(j,:)));
end;
texcode=strvcat(texcode,sprintf('\\end{array} '));
if nrow>1 | ncol>1
    texcode=strvcat(texcode,sprintf('\\right]'));
end
if ~isempty(opt.env)
    texcode=strvcat(texcode,sprintf(['\\end{',opt.env,'}']));
end

if ~isempty(opt.filename)
  eval(['fid=fopen(''',opt.filename,'.tex'',''w'')']);
  for j=1:size(texcode,1);
      rowj=texcode(j,:);
      rowj=strrep(rowj,'\','\\');
      fprintf(fid,[rowj,' \n']);
  end
  fclose(fid);
end;
