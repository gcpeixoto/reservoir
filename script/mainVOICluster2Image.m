%% mainVOICluster2Image.m - converts cluster layers to image stack
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 12nd, 2015      
%
%   Description: export cluster data to image sequence

clear all; close all;

% classes
dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEDependency;  
d.graphDataDependency;
d.VOIgraphDataDependency;

%% 

% load DRT matrix
DRT = replaceInfDRT('../mat/DRT_Field.mat');
DRT = 0*DRT;

ic = 45; jc = 68; % well
P = [14,14]; % VOI rings
drt = 14;

% load VOI structure to get cluster information
[fout,wname] = dm.createWellDir(ic,jc,'mat');
vfile = dm.setVOIFile(fout,wname,drt);
load(vfile,'VOISt');

% saving dir
svdir = '../img/';

% image format 
fmt = '.jpg';

j = 3; % component number
ind = VOISt.compVoxelInds{j};
DRT(ind) = 1;
vol = DRT( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
vol = mat2gray(vol);


for k = 1:size(vol,3);
    aux = num2str(k); 
    reldir = strcat('subvol_drt',num2str(drt),'_comp',num2str(j),'_img');
    [~,~,~] = mkdir(svdir,reldir);    
    name = strcat('subvol_drt',num2str(drt),'_comp',num2str(j),'_');    
    img = vol(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
end
    
%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
