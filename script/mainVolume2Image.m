%% mainVolume2Image.m
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 12nd, 2015      
%
%   Description: export reservoir volume data to image sequence

clear all; close all; clc;

% classes
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

% option to compress to zip (require .zip)
compr = true;

% sweep layers
for k = 1:size(PHI,3)    
    
    % ATTENTION: to save the image stack inversely:     
    %aux = num2str( size(PHI,3) - k + 1); 
    aux = num2str(k); 
    
    reldir = 'phi_img';
    [~,~,~] = mkdir(svdir,reldir); % creates dir
    name = 'spe2_phi_';    
    img = PHI(:,:,k);
        
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'kx_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'spe2_kx_';    
    img = KX(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'ky_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'spe2_ky_';    
    img = KX(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'kz_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'spe2_kz_';    
    img = KZ(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
    
    reldir = 'drt_img';
    [~,~,~] = mkdir(svdir,reldir);    
    name = 'spe2_drt_';    
    img = DRT(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
        
end

% compress files?
% Unix command: zip
if compr == true && isunix                  
        
    scal = 'phi';
    dm.zipCmd(scal,true);            
        
    scal = 'kx';
    dm.zipCmd(scal,true);               
        
    scal = 'ky';
    dm.zipCmd(scal,true);               
        
    scal = 'kz';
    dm.zipCmd(scal,true);               
        
    scal = 'drt';
    dm.zipCmd(scal,true);                       
            
end

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
