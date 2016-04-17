function [PHI,KX,KY,KZ] = buildModel(d,I,J,K)
% BUILDMODEL build the SPE Model 2 as 3D arrays.
%
%   input:
%           d: SPEDisplay class handle
%       I,J,K: SPE Model 2 grid bounds
%
%   output:
%       PHI,KX,KY,KZ: porosity, permeability 3D arrays (IxJxK)
%

% .mat files
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

%. dat files
phidat = '../dat/spe_phi.dat';
perdat = '../dat/spe_perm.dat'; 

% check existence of .dat files (MANDATORY!)
fid = fopen(phidat); assert(fid ~= -1,'Porosity file is missing!');
fid = fopen(perdat); assert(fid ~= -1,'Permeability file is missing!');

% check if code was already run
fid1 = fopen(phiname);
fid2 = fopen(kxname);
fid3 = fopen(kyname);
fid4 = fopen(kzname);
fids = [fid1,fid2,fid3,fid4]; 
status = find(fids == -1);

if isempty(status)              % reload
    d.dispMsg(d.extSPEReload);                            
    load(phiname,'PHI');    
    load(kxname,'KX');
    load(kyname,'KY');
    load(kzname,'KZ');         
else
    d.dispMsg(d.extSPERebuild); % rebuild        
    phi = load(phidat,'-ascii');
    per = load(perdat,'-ascii');
    [PHI,KX,KY,KZ] = assemble3DArrays(phi,per,I,J,K);    
        
    disp('Saving .mat files...');    
    save('../mat/PHI.mat','PHI');
    save('../mat/KX.mat','KX');
    save('../mat/KY.mat','KY');
    save('../mat/KZ.mat','KZ');
end
    
fclose('all'); 

end

