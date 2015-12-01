%% mainVOIGraphFlow

clear all; close all; clc; 

%% Load properties 

aux = load('../mat/PHI.mat');
PHI = aux.PHI;
aux = load('../mat/KN_Field.mat');
KN = aux.KN;

%% Load metrics

drtVal = '7'; % DRT to get
% metrics data structure
aux = load( strcat('../mat/VOI_DRT_',drtVal,'_MetricsData_.mat') );
metrics = aux.metrics; % components with more than 10 voxels

% DRT data structure
aux = load( strcat('../mat/DRT_VOI_',drtVal,'_Well_I45_J68.mat') );
VOISt = aux.VOISt; % all components

%for n = 1:length(metrics.idComp)
for n = 10
    
    cvc = VOISt.compVoxelCoords{n};
    cvi = VOISt.compVoxelInds{n};
                
    MadjFlowPHI = sparse( length(cvc),length(cvc) );
    MadjFlowKN = sparse( length(cvc),length(cvc) );
            
    for nv = 1:length(cvc)
        [ vc6n, vi6n, phi6n, kn6n ] = getVoxel6Neigh( cvc(nv,:), cvc, cvi, PHI, KN );        
        
        phic = PHI( cvi(nv) );
        knc = KN( cvi(nv) );
        
        for viz = 1:6
                        
            if phi6n(viz) > phic && phi6n(viz) ~= 0 
                MadjFlowPHI( vi6n(viz), nv ) = 1;
            elseif phi6n(viz) < phic && phi6n(viz) ~= 0 
                MadjFlowPHI( nv, vi6n(viz) ) = 1;
            end
            
            if kn6n(viz) > knc && kn6n(viz) ~= 0 
                MadjFlowKN( vi6n(viz), nv ) = 1;
            elseif phi6n(viz) < phic && kn6n(viz) ~= 0 
                MadjFlowKN( nv, vi6n(viz) ) = 1;
            end
            
        end
        
    end
      
    
    
    
end
