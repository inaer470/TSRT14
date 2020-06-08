function texcode=textable(V,varargin)

%Generate a LaTex table from a matlab matrix
%   texcode=textable(V,Property1,Value1,...)
%
%   V is the matrix to be tabulated.
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
%   xlabel     {''}             Array with strings of column labels
%   ylabel     {''}             Array with strings of row labels
%   title      {''}             String with table title
%
%   Examples:
%    d=1:6;
%    p=0:0.1:1;
%    for m=1:length(d)
%       xlabel{m}=['d=',num2str(d(m))];
%       for n=1:length(p)
%          ylabel{n}=['p=',num2str(p(n))];
%          A(n,m)=erfinv(chi2dist(d(m)),p(n));
%       end
%    end
%    textable(A,'ylabel',ylabel,'xlabel',xlabel,'title','h for P(chi2(d)<h)=p')
%
%   See also: texmatrix, ss/tex, tf/tex

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $


opt=struct('filename','','decimals',1,'xlabel','','ylabel','','title','');
opt=optset(opt,varargin);

[nrow,ncol]=size(V);

if isempty(opt.ylabel)
  tabarg='|';
  nylabel=0;
else
  tabarg='|l|';
  nylabel=1;
end
for i=1:ncol;
  tabarg=[tabarg,'c|'];
end;

texcode=sprintf(['\\begin{tabular}{',tabarg,'} ']);
texcode=strvcat(texcode,sprintf('\\hline'));
if ~isempty(opt.title);
  texcode=strvcat(texcode,sprintf(['\\multicolumn{',num2str(ncol+nylabel),'}{|c|}{',opt.title,'} \\\\  \\hline']));
end;

if ~isempty(opt.xlabel)
   row1=[];
   for i=1:ncol;
     row1=[row1,sprintf([' & ',opt.xlabel{i}])];
   end;
   row1=[row1,sprintf(' \\\\ \\hline \\hline')];
   texcode=strvcat(texcode,row1);
end

for j=1:nrow;
  if ~isempty(opt.ylabel)
    rowj=[opt.ylabel{j},' &  '];
  else
    rowj=[];
  end
  rowj=[rowj,'%.',num2str(opt.decimals),'f'];
  for i=2:ncol;
    rowj=[rowj,' & %.',num2str(opt.decimals),'f'];
  end;
  rowj=[rowj,' \\\\'];
  texcode=strvcat(texcode,sprintf(rowj,V(j,:)));
end;

texcode=strvcat(texcode,sprintf('\\hline  '));
texcode=strvcat(texcode,sprintf('\\end{tabular}  '));


if ~isempty(opt.filename)
  eval(['fid=fopen(''',opt.filename,'.tex'',''w'')']);
  for j=1:size(texcode,1);
      rowj=texcode(j,:);
      rowj=strrep(rowj,'\','\\');
      fprintf(fid,[rowj,' \n']);
  end
  fclose(fid);
end;
