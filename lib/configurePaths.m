%% configurePaths
% This script adds dependent folders to your Matlab path.
% However, you can also use the GUI 'Set Path'.

bdir = '/Users/gustavo/Programs/';                  % specify your base dir 
gibbon_path = fullfile(bdir,'gibbon');              % Gibbon dir
matlab_tools_path = fullfile(bdir,'matlab-tools');  % Matlab tools
snap_path = fullfile(bdir,'snap');                  % SNAP dir


% set paths
% addDependencies(path1,path2,...,pathn)    
depaths = addDependencies(gibbon_path,matlab_tools_path,snap_path);
for i = 1:length(depaths)
    
    if ~exist(depaths{i},'dir')
        warning('Path not found: %s. Check your directory layout.',depaths{i});
    else
        addpath(genpath(depaths{i}));
    end
end

% lib + script
addpath(pwd); addpath('../script/');


% clear variables
clear('i','bdir','depaths');
clear -regexp \W*path