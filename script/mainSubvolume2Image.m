%% mainSubvolume2Image.m - converts subvolume layers to image stack
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 16th, 2015      
%
%   Description: export reservoir subvolume data to image sequence

clear all; close all; clc;

dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEDependency;  
d.graphDataDependency;

%% 
% file paths
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

load(phiname,'PHI');
load(kxname,'KX');
load(kyname,'KY');
load(kzname,'KZ');

DRT = replaceInfDRT('../mat/DRT_Field.mat');

ic = 45; jc = 68; % well
P = [14,14]; % VOI rings

DRT = DRT( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
PHI = PHI( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
KX = KX( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
KY = KY( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
KZ = KZ( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );

drt = 13;
DRT = DRT == drt;

met  = 'gray';
switch met
    case 'bin'
        PHI = uint8( PHI./max( max ( max(PHI) ) ) );
        KX = uint8(KX./max( max ( max(KX) ) ) );
        KY = uint8(KY./max( max ( max(KY) ) ));
        KZ = uint8(KZ./max( max ( max(KZ) ) ));
        DRT = uint8(DRT./max( max ( max(DRT) ) ));

    case 'gray'
        PHI = mat2gray(PHI);
        KX = mat2gray(KX);
        KY = mat2gray(KY);
        KZ = mat2gray(KZ);
        DRT = mat2gray(DRT);
end

% saving dir
svdir = '../img/';

% image format 
%fmt = '.tif';
fmt = '.jpg';

% sweep layers
for k = 1:size(PHI,3)    
    
    % ATTENTION: to save the image stack inversely:     
    %aux = num2str( size(PHI,3) - k + 1); 
    aux = num2str(k); 
    
    reldir = 'subvol_phi_img';
    [~,~,~] = mkdir(svdir,reldir); % creates dir
    name = 'subvol_phi_';    
    img = PHI(:,:,k);
        
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'subvol_kx_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'subvol_kx_';    
    img = KX(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'subvol_ky_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'subvol_ky_';    
    img = KX(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'subvol_kz_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'subvol_kz_';    
    img = KZ(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = strcat('subvol_drt',num2str(drt),'_img');
    [~,~,~] = mkdir(svdir,reldir);    
    name = strcat('subvol_drt',num2str(drt),'_');    
    img = DRT(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
        
end

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
