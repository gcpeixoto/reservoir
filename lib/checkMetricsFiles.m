function matFilesOut = checkMetricsFiles(matFiles,dbase)
% CHECKMETRICSFILES checks if the metrics files was already computed. 
%                   if yes, the file structure is reduced to be rerun
%                   if no, code run normally
%
%                   matFiles: list of files (struct)
%                   dbase: dir (path to folder)

% test
assert(isstruct(matFiles),'checkMetricsFiles > argument is not a struct');

% bypasses metrics files from the list
marker = [];
for i = 3:length(matFiles)
    [~,nm,~] = fileparts(matFiles(i).name);
    
    if isempty(strfind(nm,'MetricsData')) && ...
       isempty(strfind(nm,'LinRegrData'))        
            marker = [marker; i];
    end    
end

% if any found, reduces the list
if isempty(marker)
    matFilesOut = matFiles; 
else    
    warning(strcat('Metrics files already found in',dbase,'. Rerunning...'));
    matFilesOut = matFiles(marker);
end



