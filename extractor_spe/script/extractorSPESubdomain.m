%% extractor_SPE_subdomains

%% DEFAULTS
clear all; close all; clc; format long;
diary('extractor_SPE_subdomains.log');
diary on

%% LOAD FILES

disp('Loading files...');

% files
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

fid1 = fopen(phiname);
fid2 = fopen(kxname);
fid3 = fopen(kyname);
fid4 = fopen(kzname);
if fid1 == -1 || fid2 == -1 || fid3 == -1 || fid4 == -1
    error(['Some required file was not found.'...
            'Please, run extractor again' ...
            'with option save .mat enabled']);
else
    a = load(phiname);
    PHI = a.PHI;
    
    a = load(kxname);
    KX = a.KX;
    
    a = load(kyname);
    KY = a.KY;
    
    a = load(kzname);
    KZ = a.KZ;
    
    disp('Files loaded.');  
end

ic = input('-----> Choose central voxel i coordinate: \n');
jc = input('-----> Choose central voxel j coordinate: \n');
kc = input('-----> Choose central voxel k coordinate: \n');

[PHIV, VC, VN] = getVoxelNeigh(ic,jc,kc,PHI);
[KXV, ~, ~] = getVoxelNeigh(ic,jc,kc,KX);
[KYV, ~, ~] = getVoxelNeigh(ic,jc,kc,KY);
[KZV, ~, ~] = getVoxelNeigh(ic,jc,kc,KZ);

KXVN = sqrt( KXV.^2 + KYV.^2 + KZV.^2 );

if ~all( PHIV(:) ) 
    error('Zero porosity found at the voxel neighborhood chosen. FZI will produce NaN.');    
else
    PHIVZ = PHIV./(1.0 - PHIV);
    RQIV = 0.0314*sqrt( KXVN./PHIV );
    FZIV = RQIV./PHIVZ;
    DRTV = round( 2*log( FZIV ) + 10.6 );
end

%% PLOTTING 

% 3D voxel neighbourhood 
plotVoxelNeigh3D(ic,jc,kc,DRTV,0.7);

% neighbourhood layer k = kc-1 
plotVoxelNeigh2DLayer(ic,jc,kc,DRTV,1);

% neighbourhood layer k = kc
plotVoxelNeigh2DLayer(ic,jc,kc,DRTV,2);

% neighbourhood layer k = kc+1 
plotVoxelNeigh2DLayer(ic,jc,kc,DRTV,3);

%close all
diary off