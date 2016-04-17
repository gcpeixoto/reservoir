function [PHI,KX,KY,KZ] = loadMatFiles
%LOADMATFILES Load required .mat files

aux = load('../mat/PHI.mat');
PHI = aux.PHI;
aux = load('../mat/KX.mat');
KX = aux.KX;
aux = load('../mat/KY.mat');
KY = aux.KY;
aux = load('../mat/KZ.mat');
KZ = aux.KZ;

end

