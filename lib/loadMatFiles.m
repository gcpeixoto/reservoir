function [PHI,KX,KY,KZ,KN,PHIZ,RQI,FZI,DRT] = loadMatFiles
%LOADMATFILES Load required .mat files

aux = load('../mat/PHI.mat');
PHI = aux.PHI;
aux = load('../mat/KX.mat');
KX = aux.KX;
aux = load('../mat/KY.mat');
KY = aux.KY;
aux = load('../mat/KZ.mat');
KZ = aux.KZ;
aux = load('../mat/KN.mat');
KN = aux.KN;
aux = load('../mat/PHIZ.mat');
PHIZ = aux.PHIZ;
aux = load('../mat/RQI.mat');
RQI = aux.RQI;
aux = load('../mat/FZI.mat');
FZI = aux.FZI;
aux = load('../mat/DRT.mat');
DRT = aux.DRT;

 

end

