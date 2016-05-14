%% mainVOIDRTConnections - finds clusters for the VOI based on DRT values
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27th, 2015      
%
%   Description: gets the network and connected components
%                for a VOI and stores data structure
%
% REMARK: encouraged to use non-interactive execution to not 'block' GUI.

%% Defaults

clear all; close all; clc;

% classes
dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEDependency;    

%% INPUTS

ic = 45; jc = 68; kc = 1; % seed voxel (surface)
P = [14,14,84]; % VOI rings (size) 

% voxel connectivity
% REMARK: 26-neigh is invalid for CMG (no flow; finite volume approach)
dv = setNeighDist('6'); 

%{ 
    Approaches to compute the clusters:
    The flag 'DRT_strategy' can have two values:
    i) 'reservoir': the code will compute all the connected 
                    components belonging to the DRT values found
                    only over the central well of the reservoir.
   ii) 'well':      ditto, except that the computation will consider
                    the DRT values over all the reservoir.

%}
DRT_strategy = 'reservoir';
%DRT_strategy = 'well';

%% LOAD

% load DRT
[~,~,~,~,~,~,~,~,DRT] = loadMatFiles;

% get VOI
[drt,VN,IND] = getVoxelRegion(ic,jc,kc,P,DRT);
nvres = length(drt); % total number of voxels in the VOI

switch DRT_strategy
    
    case 'well'    
        
        drtVec = DRT(ic,jc,:); 
        drtVec = unique ( reshape(drtVec, [size(DRT,3), 1]) );
        
        % creates output mat dir 
        [dirp,wname] = dm.createWellDir(ic,jc,'mat','w');
        
    case 'reservoir'
        
        drtVec = unique(drt(drt>0)); % bypass DRT = 0        
        
        % creates output mat dir 
        [dirp,wname] = dm.createWellDir(ic,jc,'mat','r');
end

% defining structure fields
drtVOI.seed = [ic,jc,kc];
drtVOI.rings = P;

% sweeping DRT values of the well
for i = 1:length(drtVec)
    
    drtval = drtVec(i); % DRT
        
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
        
    fprintf('----> Computing DRT = %d... \n',drtVOI.value{m});
    
    coordsDRT = drtVOI.voxels{m};
    indz = drtVOI.voxelNodes{m};  
    
    disp('----> Computing distances...');        
    indIJ = [];
    for i = 1:size(coordsDRT,1)        
        for j = i:size(coordsDRT,1) % d(i,j) = d(j,i)                                                            
              if i ~= j % skipping null distance 
                  dist = sqrt( ( coordsDRT(i,1) - coordsDRT(j,1) )^2 + ...
                               ( coordsDRT(i,2) - coordsDRT(j,2) )^2 + ...
                               ( coordsDRT(i,3) - coordsDRT(j,3) )^2 ); 
                  
                  % detecting neighbour voxels                  
                  if dist <= dv     % connectivity criterion
                      indIJ = [ indIJ; [ i j ] ];
                      edgeList = indIJ;                       
                  end                  
              end
         end
    end   
    
    if isempty(indIJ)
        fprintf('----> No connections found for DRT = %d... \n',drtVOI.value{m});
        continue;
    end
    
    aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
    indIJ = [ indIJ; aux ];          % filling
 
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
    VOISt.allNVoxels = nvres;                % total number of voxels in the reservoir
    
    for idcomp = 1:ncomp
        
        fprintf('--------> Computing component = %d... \n',idcomp);
        
        aux = coordsDRT( members{idcomp},: );
        VOISt.compVoxelCoords{idcomp} = aux; 
        VOISt.compVoxelInds{idcomp} = indz( members{idcomp} );        
        VOISt.compNNodes{idcomp} = compSizes(idcomp);          
        
    end
    
    % save .mat file
    fname = strcat(dirp,'VOI_DRT_',num2str( drtVOI.value{m} ),'_',...
                          wname,'.mat');
    save(fname,'VOISt');
    fprintf('----> %s saved. \n',fname);
    
    clear VOISt
end

%% ENDINGS
d.printings(d.progStat{2});
dm.deactivateLog;
