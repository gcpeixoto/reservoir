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
kzname  = '../mat/KZ.mat';

fid1 = fopen(phiname);
fid2 = fopen(kxname);
fid3 = fopen(kzname);
if fopen(fid1) == -1
    error('Porosity file not found.');
end
if fopen(fid2) == -1
    error('Permeability - kx file not found.');
end
if fopen(fid3) == -1
    error('Permeability - kz file not found.');
end

phi = load(phiname);
kx = load(kxname);
kz = load(kzname);

disp('Files loaded.');  

ri = 1;
rj = 2;
rk = 3;

iref = 10;
jref = 20;
kref = 30;

I = 60; J = 220; K = 85;

if ~( ri <= I - iref ) || ~( rj <= J - jref ) || ~( rk <= K - kref )
    error('Subdomain radius exceeded in some direction.');
end

% Reference point
kx.KX(iref,jref,kref);


diary off