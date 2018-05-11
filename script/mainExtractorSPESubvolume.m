%% mainExtractorSPESubvolume.m - Subvolume data extractor for oil wells
%      
%     authors: Dr. Gustavo Peixoto de Oliveira
%              Dr. Waldir Leite Roque
%              @Federal University of Paraiba
%     mail: gustavo.oliveira@ci.ufpb.br    
%     date: Sep 21st, 2015  
%
% Description: extracts a subvolume of the petroleum field for 
%              local analysis of quantities
%              
%              User enters with: 
%              (ic,jc,kc), the seed voxel coordinates and
%                       P, the radius of the Moore's neighbourhood
%              to get a subvolume of (2*P+1)^3 voxels

%% DEFAULTS

clear all; close all; clc;

% classes
dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEDependency;    

%% LOAD FILES

% file paths
[PHI,KX,KY,KZ,~,~,~,~,~] = loadMatFiles;

%% INPUT DATA 

% plot selection
plt_drt = false;              % flag to plot DRTs in the subvolume 

ic = input(d.dispCCoord(1));
jc = input(d.dispCCoord(2));
kc = input(d.dispCCoord(3));
P  = input(d.extSPESP);

%% SUBVOLUME

[PHIV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,PHI);
[KXV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KX);
[KYV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KY);
[KZV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KZ);

KNV = sqrt( KXV.^2 + KYV.^2 + KZV.^2 );     % permeability norm
PHIVZ = PHIV./(1.0 - PHIV);                 % normalized porosity
RQIV = 0.0314*sqrt( KNV./PHIV );            % RQI
FZIV = RQIV./PHIVZ;                         % FZI
DRTV = round( 2*log( FZIV ) + 10.6 );       % DRT

%% PLOTTING 

% 3D voxel neighbourhood 
plotVoxelNeigh3D(ic,jc,kc,P,PHIVZ,1.0);

% plot all DRTs in the subvolume
if plt_drt == true
    drt = unique(DRTV(:));
    for i = 1:length(drt)
        plotVoxelNeighByValue(ic,jc,kc,P,DRTV,drt(i),1.0,'DRT',false);   
    end 
end
  
%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;