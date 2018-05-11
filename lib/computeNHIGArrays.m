function [ RQI_NHIG, FZI_NHIG, DRT_NHIG ] = computeNHIGArrays(PHI,KN,SWIR,m)
% COMPUTENHIGARRAYS compute 3D arrays for RQI, FZI and DRT 
%                   according to the NHIG model for the SPE model
%
%   input:
%          PHI: porosity
%           KN: SPE norm-permeability 3D array
%         SWIR: irreducible water saturation 3D array
%            m: cementation factor (scalar)
%

assert(m>=1,'m must be >=1'); % cementation factor

% NH-IG model (Tati's model)
 

kIG = KN./((1-SWIR).^3);    % permeability model of IG
PHIZ = PHI./(1-PHI);    
RQI_NHIG = 0.0314*sqrt( kIG./PHI );
FZI_NHIG = RQI_NHIG./( PHIZ.*PHI.^(m-1) );

[~,coords,~] = maskVoxels( FZI_NHIG,Inf );
    if ~isempty( coords ) 
        for e = 1:size(coords,1)
            FZI_NHIG( coords(e,1), coords(e,2), coords(e,3) ) = 0.0;
        end
    end

DRT_NHIG = round( 2*log( FZI_NHIG ) + 10.6 );        
DRT_NHIG = replaceInfDRT(DRT_NHIG); % replaces DRT = Inf with DRT = 0

% saving 3D arrays                
%disp('Saving .mat files...');    
    
%save('../mat/RQI_NHIG.mat','RQI_NHIG'); 
%save('../mat/FZI_NHIG.mat','FZI_NHIG'); 
%save('../mat/DRT_NHIG.mat','DRT_NHIG'); 


end

