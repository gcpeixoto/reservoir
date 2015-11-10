%% mainVOIConnections
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27st, 2015      
%
%   Description: gets the network and connected components
%                for a VOI around a given seed voxel 
%

clear all; close all; clc;

%%

% load DRT matrix
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;

ic = 45; jc = 68; kc = 1; % seed voxel
P = [14,14,84]; % VOI rings
 
% DRT values over the well
drtWell = DRT( ic, jc,:); 
drtWell = unique ( reshape(drtWell, [size(DRT,3) 1]) );

% get VOI
[drt,VN,IND] = getVoxelRegion(ic,jc,kc,P,DRT);

%%%% PLOT VOI
% [F,V,C]=ind2patch(IND,DRT,'v');
% hold on
% patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.5);
% axis equal; view(3); axis tight; axis vis3d; grid off;
% set(gca,'ZDir','reverse');
% print('-dpdf','-r0','../figs/graphpath/VOI');
%%%%

% defining structure fields
drtVOI.seed = [ic,jc,kc];
drtVOI.rings = P;

% sweeping DRT values of the well
for i = 1:length(drtWell)
    
    drtval = drtWell(i); % DRT
    
    I = find( drt == drtval ); % getting local indices
    v = [];
    ind = [];
    for j = 1:length(I)
        v = [ v; VN(I(j),:) ]; % getting global voxels
        ind = [ ind; sub2ind( size(DRT),VN(I(j),1),VN(I(j),2),VN(I(j),3) ) ]; % getting global indices
    end
    
    % storing in the structure
    drtVOI.value{i} = drtval;
    drtVOI.voxels{i} = v;
    drtVOI.voxelNodes{i} = ind;       
    
end


% loop to get VOI graph data 
for m = 1:length( drtVOI.value )
    
    fprintf('----> m = %d... \n',m); 
    
    coordsDRT = drtVOI.voxels{m};
    indz = drtVOI.voxelNodes{m};  
    
    disp('----> Computing distances...');        
    indIJ = [];
    for i = 1:size(coordsDRT,1)
        fprintf('----> i = %d... \n',size(coordsDRT,1)-i); 
        for j = i:size(coordsDRT,1)
                                                            
              if i ~= j % skipping null distance 
                  dist = sqrt( ( coordsDRT(i,1) - coordsDRT(j,1) )^2 + ...
                               ( coordsDRT(i,2) - coordsDRT(j,2) )^2 + ...
                               ( coordsDRT(i,3) - coordsDRT(j,3) )^2 ); 
                  
                  % detecting neighbour voxels
                  if dist <= sqrt(3)
                      indIJ = [ indIJ; [ i j ] ];
                      edgeList = indIJ;                       
                  end                  
              end
         end
    end   
    aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
    indIJ = [ indIJ; aux ]; % filling
 
    disp('----> Computing adjacency matrix...');            
    % creates adjacency matrix n x n by marking 1 for connected nodes
    MDadj = sparse( indIJ(:,1),indIJ(:,2),1,size(coordsDRT,1),size(coordsDRT,1) ); 
            
    disp('----> Finding network components...');
    [ncomp,compSizes,members] = networkComponents(MDadj);
    
    % storing in strucuture
    VOISt.value = drtVOI.value{m};           % DRT value
    VOISt.allAdjMatrix = MDadj;              % graph adjacency matrix
    VOISt.allAdjEdgeList = edgeList;         % graph edge list
    VOISt.allVoxelCoords = coordsDRT;        % voxel coordinates (i,j,k)     
    VOISt.allVoxelInds = indz;               % voxel linear indices
    VOISt.allNComps = ncomp;                 % number of components in the graph
    
    for idcomp = 1:ncomp
        
        fprintf('--------> Computing component = %d... \n',idcomp);
        
        aux = coordsDRT( members{idcomp},: );
        VOISt.compVoxelCoords{idcomp} = aux; 
        VOISt.compVoxelInds{idcomp} = indz( members{idcomp} );        
        VOISt.compNNodes{idcomp} = compSizes(idcomp);          
        
    end    
    
    save( strcat('../mat/DRT_VOI_',num2str( drtVOI.value{m} ),...
                          '_Well_I',num2str( ic ),... 
                              '_J',num2str( jc ),'.mat'),'VOISt'); % saving
    disp('----> .mat file saved.')
end





