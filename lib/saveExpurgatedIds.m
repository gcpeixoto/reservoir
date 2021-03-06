function saveExpurgatedIds(idvec,i,j,wdir)
% SAVEEXPURGATEDIDS write to file the well's 
%      depth z-coordinates where a value of 0 was
%      found for porosity. These values are removed in the analysis
%      and this file is only intended to identify such points.
%   
%       input: 
%              idvec: z coordinates array (nz x 1)
%                i,j: well coordinates
%               wdir: well dir

aux = [];
for k=1:length(idvec)    
    aux = [ aux; i j idvec(k) ];
    fname = strcat(wdir,'ExpurgatedPoints_I_',num2str(i),'J_',num2str(j),'.csv');    
    dlmwrite(fname,aux);
        
end

