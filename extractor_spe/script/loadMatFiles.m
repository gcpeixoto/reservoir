function [ PHI,KX,KY,KZ ] = loadMatFiles( phiname,kxname,kyname,kzname )
%   LOADMATFILES loads .mat files
%

disp('Loading files...');

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

end

