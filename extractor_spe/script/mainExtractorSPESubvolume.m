%% mainExtractorSPESubvolume.m - Subvolume data extractor for oil wells
      
%     authors: Dr. Gustavo Peixoto de Oliveira
%              Dr. Waldir Leite Roque
%              @Federal University of Paraiba
%     mail: gustavo.oliveira@ci.ufpb.br    
%     date: Sep 21st, 2015  

%% DEFAULTS

clear all; close all; clc;
activateLog(mfilename);
setOptions;
splshScreenSub;

%% LOAD FILES

% file paths
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

load(phiname,'PHI');
load(kxname,'KX');
load(kyname,'KY');
load(kzname,'KZ');

%% INPUT DATA 

ic = input('-----> Choose central voxel i coordinate: \n');
jc = input('-----> Choose central voxel j coordinate: \n');
kc = input('-----> Choose central voxel k coordinate: \n');
P  = input('-----> Choose P ring radius: \n');

%% COMPUTATION

[PHIV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,PHI);
[KXV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KX);
[KYV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KY);
[KZV,~,~,~] = getVoxelNeighRing(ic,jc,kc,P,KZ);

KXVN = sqrt( KXV.^2 + KYV.^2 + KZV.^2 );
PHIVZ = PHIV./(1.0 - PHIV);
RQIV = 0.0314*sqrt( KXVN./PHIV );
FZIV = RQIV./PHIVZ;
DRTV = round( 2*log( FZIV ) + 10.6 );

%% PLOTTING 

% 3D voxel neighbourhood 
plotVoxelNeigh3D(ic,jc,kc,P,PHIVZ,1.0);

% plot all DRTs
drt = unique(DRTV(:));
for i = 1:length(drt)
     plotVoxelNeighByValue(ic,jc,kc,P,DRTV,drt(i),1.0,'DRT');   
end    
  
diary off