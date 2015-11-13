%% mainPlotPipeline
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27th, 2015      
%
%   Description: gets the network and connected components
%                of the whole field for a given DRT and plots 
%                a pipeline mimic to feature the oil duct
%                path from the voxel C to F, where 
%                C is the vertex of maximum closeness centrality and 
%                F is the vertex of minimum closeness centrality

clear all; close all; clc

%% Load data

% DRT
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

drtVal = '6'; % DRT to get

% metrics data structure
aux = load( strcat('../mat/DRT_',drtVal,'_MetricsData_.mat') );
metrics = aux.metrics;

% DRT data structure
aux = load( strcat('../mat/DRT_',drtVal,'.mat') );
drtSt = aux.drtSt;
val = drtSt.value;

% data of components with high performance of correlation from metrics data
ncomp = numel(metrics.idComp);      % number of components
compVoxelInds = cell(ncomp,1);      % global linear indices of component's voxels
compVoxelCoords = cell(ncomp,1);    % global voxel coords of component's voxels
adjMat = cell(ncomp,1);             % component's adjacency matrix
closeness = cell(ncomp,1);          % closeness centrality
idclo = cell(ncomp,1);              % min closeness voxel indices
idfar = cell(ncomp,1);              % max farness voxel indices
distsToCloseness = cell(ncomp,1);   % cell of distances to the closest point
pathsToCloseness = cell(ncomp,1);   % cell of path local indices to the closest point
pathsVoxelCoords = cell(ncomp,1);   % cell of path global voxels to the closest point

% --------- component loop

for n = 1:ncomp                     
    
    compVoxelInds{n} = drtSt.compVoxelInds{ metrics.idComp{n} };        
    compVoxelCoords{n} = drtSt.compVoxelCoords{ metrics.idComp{n} };    
    adjMat{n} = metrics.adjMatrix{n};
    closeness{n} = metrics.closenessCentrality{n};
    idclo{n} = find( closeness{n} == max( closeness{n} ) ); % most central node local id
    idfar{n} = find( closeness{n} == min( closeness{n} ) ); % most distant node local id
          
    % More than one voxel might have maximum(minimum) closeness. These 'ifs' 
    % choose the first found in the list 
    if numel(idclo{n}) > 1, idclo{n} = idclo{n}(1); end
    if numel(idfar{n}) > 1, idfar{n} = idfar{n}(1); end
    
    [dists,paths] = dijkstra( adjMat{n},idfar{n},idclo{n} ); % path from the closest to more distant node in the network
    distsToCloseness{n} = dists;
    pathsToCloseness{n} = paths;  
    
    inds{n} = compVoxelInds{n}( pathsToCloseness{n} );       % global lin. indices of the paths
    coords{n} = compVoxelCoords{n}( pathsToCloseness{n},: ); % global voxel coords of the paths
    
    %----- PLOT PIPELINE
    % view([0,90]) - projection
        
    plotVoxelGraphComp(DRT,compVoxelInds{n},val,n,0.03,'gray','no', [-46,28] ); % plot whole connected comp    
    %plotVoxelGraphComp(DRT,inds{n},val,n,0.03,'gray','no', [-46,28] ); % plot only the voxels of the path
    
    % plots the oriented path
    hold on
    h = plot3(coords{n}(:,2), coords{n}(:,1), coords{n}(:,3));
    set( h,'Color',rgb('FireBrick'),'LineWidth',3,'Marker','o',...
      'MarkerSize',10,'MarkerFaceColor',rgb('ForestGreen'),...
      'MarkerEdgeColor',rgb('ForestGreen'));
    
  % marks the closest point
    hold on
    h2 = plot3( compVoxelCoords{n}( idclo{n},2 ),compVoxelCoords{n}( idclo{n},1 ),compVoxelCoords{n}( idclo{n},3 ) );
    set( h2,'Marker','o',...
      'MarkerSize',25,'MarkerFaceColor',rgb('Gold'),...
      'MarkerEdgeColor',rgb('Gold'));    
  
    fname = strcat('../figs/graphPath/pipeline_DRT_',num2str(drtVal),'_Comp_',num2str(n) );
    %print('-depsc2','-r0',fname);
    print('-dpdf','-r0',fname);
        
end
