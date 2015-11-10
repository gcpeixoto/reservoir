%% mainPlotPipeline

clear all; close all; clc

%% Load data

% DRT
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

% metrics
aux = load('../mat/DRT_6_MetricsData_.mat');
metrics = aux.metrics;

% structure
aux = load('../mat/DRT_6.mat');
drtSt = aux.drtSt;
val = drtSt.value;

% getting high performance components from metrics data
ncomp = numel(metrics.idComp);
compVoxelInds = cell(ncomp,1);
compVoxelCoords = cell(ncomp,1);
adjMat = cell(ncomp,1);
closeness = cell(ncomp,1);
idclo = cell(ncomp,1); % min closeness voxel indices
idfar = cell(ncomp,1); % max farness voxel indices
distsToCloseness = cell(ncomp,1);
pathsToCloseness = cell(ncomp,1);
pathsVoxelCoords = cell(ncomp,1);

for n = 1:ncomp
    compVoxelInds{n} = drtSt.compVoxelInds{ metrics.idComp{n} };
    compVoxelCoords{n} = drtSt.compVoxelCoords{ metrics.idComp{n} };
    adjMat{n} = metrics.adjMatrix{n};
    closeness{n} = metrics.closenessCentrality{n};
    idclo{n} = find( closeness{n} == max( closeness{n} ) );
    idfar{n} = find( closeness{n} == min( closeness{n} ) );
        
    if numel(idclo{n}) > 1
        idclo{n} = idclo{n}(1); % gets only 1
    end
    if numel(idfar{n}) > 1
        idfar{n} = idfar{n}(1);
    end
    
    [dists,paths] = dijkstra( adjMat{n},idfar{n},idclo{n} );
    distsToCloseness{n} = dists;
    pathsToCloseness{n} = paths;  
    
    inds{n} = compVoxelInds{n}( pathsToCloseness{n} );
    coords{n} = compVoxelCoords{n}( pathsToCloseness{n},: );
    
    plotVoxelGraphComp(DRT,compVoxelInds{n},val,n,0.03,'gray','no', [-46,28] );
    hold on
    %plotVoxelGraphComp(DRT,inds{n},val,n,0.03,'gray','no', [-46,28] );
    hold on
    h = plot3(coords{n}(:,2), coords{n}(:,1), coords{n}(:,3));
    set( h,'Color',rgb('FireBrick'),'LineWidth',3,'Marker','o',...
      'MarkerSize',10,'MarkerFaceColor',rgb('ForestGreen'),...
      'MarkerEdgeColor',rgb('ForestGreen'));
    hold on
    h2 = plot3( compVoxelCoords{n}( idclo{n},2 ),compVoxelCoords{n}( idclo{n},1 ),compVoxelCoords{n}( idclo{n},3 ) );
    set( h2,'Marker','o',...
      'MarkerSize',25,'MarkerFaceColor',rgb('Gold'),...
      'MarkerEdgeColor',rgb('Gold'));    
        
end


% plotVoxelGraphComp(DRT,vxinds1,6,1,0.03,'gray','no', [0,90] );
% hold on
% %scatter3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3),200,'filled');
% h=plot3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3));
% set( h,'Color',rgb('FireBrick'),'LineWidth',3,'Marker','o',...
%       'MarkerSize',10,'MarkerFaceColor',rgb('ForestGreen'),...
%       'MarkerEdgeColor',rgb('ForestGreen'));
% hold on
% h2 = plot3( vxcord1( clo1(1),2 ),vxcord1( clo1(1),1 ),vxcord1( clo1(1),3 ) );
% set( h2,'Marker','o',...
%       'MarkerSize',25,'MarkerFaceColor',rgb('Gold'),...
%       'MarkerEdgeColor',rgb('Gold'));
% 
% print('-depsc2','-r0','../figs/graphPath/pipeline1');
%   
% plotVoxelGraphComp(DRT,vxinds1,6,1,0.03,'gray','no', [-46,28] );
% hold on
% %scatter3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3),200,'filled');
% h=plot3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3));
% set( h,'Color',rgb('FireBrick'),'LineWidth',3,'Marker','o',...
%       'MarkerSize',10,'MarkerFaceColor',rgb('ForestGreen'),...
%       'MarkerEdgeColor',rgb('ForestGreen'));
% hold on
% h2 = plot3( vxcord1( clo1(1),2 ),vxcord1( clo1(1),1 ),vxcord1( clo1(1),3 ) );
% set( h2,'Marker','o',...
%       'MarkerSize',25,'MarkerFaceColor',rgb('Gold'),...
%       'MarkerEdgeColor',rgb('Gold'));
%   
% print('-depsc2','-r0','../figs/graphPath/pipeline2');
%   
% plotVoxelGraphComp(DRT,vxinds1,6,1,0.03,'gray','no', [57,54] );
% hold on
% %scatter3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3),200,'filled');
% h=plot3(pathcoords(:,2), pathcoords(:,1), pathcoords(:,3));
% set( h,'Color',rgb('FireBrick'),'LineWidth',3,'Marker','o',...
%       'MarkerSize',10,'MarkerFaceColor',rgb('ForestGreen'),...
%       'MarkerEdgeColor',rgb('ForestGreen'));
% hold on
% h2 = plot3( vxcord1( clo1(1),2 ),vxcord1( clo1(1),1 ),vxcord1( clo1(1),3 ) );
% set( h2,'Marker','o',...
%       'MarkerSize',25,'MarkerFaceColor',rgb('Gold'),...
%       'MarkerEdgeColor',rgb('Gold'));  
%   
% print('-depsc2','-r0','../figs/graphPath/pipeline3');
