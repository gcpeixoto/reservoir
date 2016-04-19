%% configurePaths
% This script adds dependent folders to your Matlab path.
% However, you can also use the GUI 'Set Path'.

resdir = '/Users/gustavo/reservoir/';                % reservoir dir
%resdir = '/home/tatiana/reservoir/';

bdir = '/Users/gustavo/Programs/';                  % third-party base dir 
%bdir = '/home/tatiana/Programs/';

gibbon_path = fullfile(bdir,'gibbon');              % Gibbon dir
matlab_tools_path = fullfile(bdir,'matlab-tools');  % Matlab tools
snap_path = fullfile(bdir,'snap');                  % SNAP dir

script_dir = fullfile(resdir,'script');             % script dir
lib_dir = fullfile(resdir,'lib');                   % lib dir

% set paths
% addDependencies(path1,path2,...,pathn)    
depaths = addDependencies(gibbon_path,matlab_tools_path,snap_path,...
                          script_dir,lib_dir);
for i = 1:length(depaths)
    
    if ~exist(depaths{i},'dir')
        warning('Path not found: %s. Check your directory layout.',depaths{i});
    else
        addpath(genpath(depaths{i}));
    end
end

% clear variables
clear('i','bdir','resdir','depaths');
clear -regexp \W*path