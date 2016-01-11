%% mainVOIMetricsAnalyzer


clear all; close all; clc; format long;
set(0,'DefaultAxesFontSize',18);  
%% 

% DRT
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

aux = load('../mat/PHI.mat');
PHI = aux.PHI;

wellfile = '_Well_I45_J68';

drtVal = '7'; % DRT to get

% metrics data structure
aux = load( strcat('../mat/VOI_DRT_',drtVal,'_MetricsData_.mat') );
metrics = aux.metrics;

% DRT data structure
aux = load( strcat('../mat/DRT_VOI_',drtVal,wellfile,'.mat') );
VOISt = aux.VOISt;
val = VOISt.value;

ncomp = numel(metrics.idComp);      % number of components

nc = 3;
%nc = ncomp;
sumcomp = 0; % total number of voxels for the current DRT (nc components)
for n = 1:nc
    
    compVoxelInds{n} = VOISt.compVoxelInds{ metrics.idComp{n} };        
    compVoxelCoords{n} = VOISt.compVoxelCoords{ metrics.idComp{n} };    
    adjMat{n} = metrics.adjMatrix{n};
    
    % centralities
    degree{n} = metrics.degreeCentrality{n};
    betweeness{n} = metrics.betweenessCentrality{n};
    closeness{n} = metrics.closenessCentrality{n};
    
    % max deg
    maxdeg{n} = max( degree{n} );                    
    imaxdeg{n} = find( degree{n} == maxdeg{n} );     
    cvcmaxdeg{n} = compVoxelCoords{n}(imaxdeg{n},:); 
    
    % min deg
    mindeg{n} = min( degree{n} );    
    imindeg{n} = find( degree{n} == mindeg{n} );    
    cvcmindeg{n} = compVoxelCoords{n}(imindeg{n},:);
    
    % max bet
    maxbet{n} = max( betweeness{n} );
    imaxbet{n} = find( betweeness{n} == maxbet{n} );
    cvcmaxbet{n} = compVoxelCoords{n}(imaxbet{n},:);
    
    % min bet
    minbet{n} = min( betweeness{n} );    
    iminbet{n} = find( betweeness{n} == minbet{n} );    
    cvcminbet{n} = compVoxelCoords{n}(iminbet{n},:);
    
    % max clo
    maxclo{n} = max( closeness{n} );
    imaxclo{n} = find( closeness{n} == maxclo{n} );    
    cvcmaxclo{n} = compVoxelCoords{n}(imaxclo{n},:);
    degmaxclo{n} = degree{n}(imaxclo{n});
    
    % min clo
    minclo{n} = min( closeness{n} );
    iminclo{n} = find( closeness{n} == minclo{n} );
    cvcminclo{n} = compVoxelCoords{n}(iminclo{n},:);
    degminclo{n} = degree{n}(iminclo{n});
                
    sumcomp = sumcomp + VOISt.compNNodes{n};
    
    % farthest voxels
    [D2CFar,CVCFar,ilims,jlims,klims] = getDists2Point(compVoxelCoords{n},cvcmaxclo,1);    
    compD2CFar{n} = D2CFar;
    compCVCFar{n} = CVCFar;
    compILims{n} = ilims;
    compJLims{n} = jlims;
    compKLims{n} = klims;
       
    plotMetricField(compVoxelCoords{n},closeness{n},'clo',n,val);
    hold on
    % highlight farthest 
    scatter3(CVCFar(:,2),CVCFar(:,1),CVCFar(:,3),300,'k');
    
    %plotMetricField(compVoxelCoords{n},betweeness{n},'bet',n,val);
    %plotMetricField(compVoxelCoords{n},degree{n},'deg',n,val);    
    
    
    
            
end
