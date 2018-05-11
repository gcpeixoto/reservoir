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
    
    KN = sqrt( KX.^2 + KY.^2 + KZ.^2 );
    PHIZ = PHI./(1.0 - PHI);
    RQI = 0.0314*sqrt( KN./PHI );
    FZI = RQI./PHIZ;

    [~,coords,~] = maskVoxels( FZI,Inf );
    if ~isempty( coords ) 
        for e = 1:size(coords,1)
            FZI( coords(e,1), coords(e,2), coords(e,3) ) = 0.0;
        end
    end

    DRT = round( 2*log( FZI ) + 10.6 );        
    DRT = replaceInfDRT(DRT); % replaces DRT = Inf with DRT = 0

    % saving 3D arrays                
    disp('Saving .mat files...');    
    
    save('../mat/PHI.mat','PHI');
    save('../mat/KX.mat','KX');
    save('../mat/KY.mat','KY');
    save('../mat/KZ.mat','KZ');    
    save('../mat/KN.mat','KN'); 
    save('../mat/PHIZ.mat','PHIZ'); 
    save('../mat/RQI.mat','RQI'); 
    save('../mat/FZI.mat','FZI'); 
    save('../mat/DRT.mat','DRT'); 
    
end
    
fclose('all'); 

end

