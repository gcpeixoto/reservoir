function DRTFiles = getDRTFile(matFiles,dbase,drtVal)
% GETDRTFILE get files from specified DRT 
%
%                   matFiles: list of files (struct)
%                   dbase: matFiles dir (char)
%                   drtVal: DRT value (int)
%                   spec: 'Well','MetricsData','' 

% test
assert(isstruct(matFiles),'getDRTFile > argument is not a struct');

for i = 3:length(matFiles)
    [~,nm,~] = fileparts(matFiles(i).name);    
    if ~isempty(strfind(nm,strcat('DRT_',num2str(drtVal),'_Well')))       
            DRTFiles{1} = strcat(dbase,matFiles(i).name);
    end
    
    if ~isempty(strfind(nm,strcat('DRT_',num2str(drtVal),'_MetricsData')))
            DRTFiles{2} = strcat(dbase,matFiles(i).name);
    end
    
    if ~isempty(strfind(nm,strcat('DRT_',num2str(drtVal),'_LinRegr')))
            DRTFiles{3} = strcat(dbase,matFiles(i).name);
    end
    
end



