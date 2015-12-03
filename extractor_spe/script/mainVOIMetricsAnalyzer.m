%% mainVOIMetricsAnalyzer


clear all; close all; clc; format long;
set(0,'DefaultAxesFontSize',18);  
%% 

% DRT
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

drtVal = '13'; % DRT to get

% metrics data structure
aux = load( strcat('../mat/VOI_DRT_',drtVal,'_MetricsData_.mat') );
metrics = aux.metrics;

% DRT data structure
aux = load( strcat('../mat/DRT_VOI_',drtVal,'_Well_I45_J68.mat') );
VOISt = aux.VOISt;
val = VOISt.value;

ncomp = numel(metrics.idComp);      % number of components

for n = 1:3%ncomp    
    
    compVoxelInds{n} = VOISt.compVoxelInds{ metrics.idComp{n} };        
    compVoxelCoords{n} = VOISt.compVoxelCoords{ metrics.idComp{n} };    
    adjMat{n} = metrics.adjMatrix{n};
    
    degree{n} = metrics.degreeCentrality{n};
    betweeness{n} = metrics.betweenessCentrality{n};
    closeness{n} = metrics.closenessCentrality{n};
    
    maxdeg{n}  = max( degree{n} );
    imaxdeg{n} = find( degree{n} == maxdeg{n} );
    cvcmaxdeg{n} = compVoxelCoords{n}(imaxdeg{n},:);
    
    maxbet{n} = max( betweeness{n} );
    imaxbet{n} = find( betweeness{n} == maxbet{n} );
    cvcmaxbet{n} = compVoxelCoords{n}(imaxbet{n},:);
    
    maxclo{n} = max( closeness{n} );
    imaxclo{n} = find( closeness{n} == maxclo{n} );
    cvcmaxclo{n} = compVoxelCoords{n}(imaxclo{n},:);
    
    plotMetricField(compVoxelCoords{n},closeness{n},'clo',n,val);
    plotMetricField(compVoxelCoords{n},betweeness{n},'bet',n,val);
    plotMetricField(compVoxelCoords{n},degree{n},'deg',n,val);
    
end