function optsigsys
%OPTSIGSYS creates a global variable with default values

global optsigsys
optsigsys.fields={'linewidth','fontsize','markersize'};
optsigsys.defaults={'1','12','5'};
for k=1:length(optsigsys.fields)
	eval(['optsigsys.',optsigsys.fields{k},'=',optsigsys.defaults{k},';'])
end
