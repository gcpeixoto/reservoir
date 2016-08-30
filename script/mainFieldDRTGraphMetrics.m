%% mainFieldDRTGraphMetrics - cluster metrics for the whole field
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27th, 2015        
%             
%   description: compute metrics and regression analysis for the whole 
%                field's networks based on the DRT field value.
%
%   requirements:
%        - pre-computed .mat files
%        - Matlab third-party additional functions
%

%% DEFAULTS
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

nofn = 10;   % minimum number of voxels to consider per component 
seps = 0.05; % linear regression epsilon for slope [1-seps,1+seps]
R2min = 0.9; % minimum R2 coefficient acceptable

%% LOAD FILES 

% load DRT
[~,~,~,~,~,~,~,~,DRT] = loadMatFiles;

matFiles = dir('../mat/DRT_*.mat'); 
numfiles = length(matFiles);
    
% sweeping DRTs
for k = 1:numfiles 
    
    load( strcat('../mat/',matFiles(k).name),'drtSt'); 
    val = drtSt.value; 
    
    fprintf('----> Sweeping DRT: %d... \n',val);
    
    avc = drtSt.allVoxelCoords;
    ncomps = drtSt.allNComps;    
    Madj = drtSt.allAdjMatrix;    
    
    metrics.drtValue = val;
    linregr.drtValue = val;
    
    count = 0;
    for idComp = 1:ncomps        
        
        cnn = drtSt.compNNodes{idComp};        
                        
        if cnn >= nofn % significative components 
            
            cvc = drtSt.compVoxelCoords{idComp};
            cvi = drtSt.compVoxelInds{idComp};            
  
            % performs linear regression
            logPHIZ = drtSt.compLogPHIZ{idComp};
            logRQI  = drtSt.compLogRQI{idComp};
            [ R, m, b ] = regression( logPHIZ, logRQI, 'one' );
                                                                            
            count = count + 1; % component counter

            %------------------ subgraph (connected component network)
            % finding vertices in the big adjacency matrix to set up the 
            % adjacency matrix for the connected component and, then, 
            % set the subgraph of the network                
            v = [];
            for e = 1:size(cvc,1)
                id = strmatch( cvc(e,:), avc );       %#ok<*MATCH2>
                v = [ v; id ];                        % global indices
            end                    
            MadjComp = subgraph( Madj, v );           % component's adjacency matrix                

            %------------------ centrality metrics 
            fprintf('----> Computing metrics for %d nodes... \n', size(v,1) );

            % SNAP interface
            edfile = saveAdjEdges(MadjComp);  
            ! ./../cpp/graphMetrics
            [nodeID,deg,clns,betw] = getMetricsData(edfile);                

            maxC = max(clns);                         % max closeness = min farness
            iC = find( clns == maxC );                % network closer nodes
            iCnode = nodeID(iC);                      % getting node id (not always == iC)
            ivC = avc( v(iCnode),: );                 % global voxel coordinates

            disp('----> Storing structures...');
            % store good components         
            metrics.idComp{count} = idComp;                
            metrics.degreeCentrality{count} = deg;
            metrics.closenessCentrality{count} = clns;
            metrics.betweenessCentrality{count} = betw;                
            metrics.maxClosenessVoxelCoords{count} = ivC;
            metrics.adjMatrix{count} = MadjComp;

            linregr.idComp{count} = idComp;
            linregr.Pearson{count} = R*R;
            linregr.slope{count} = m;
            linregr.offset{count} = b;
            linregr.logPHIZ{count} = logPHIZ;
            linregr.logRQI{count} = logRQI;
            
            % linear regression criteria (performance)                 
            if ( m >= 1-seps && m <= 1+seps ) && (R*R >= R2min)
                linregr.performance{count} = 1; % high-performance 
            else
                linregr.performance{count} = 0; % low-performance
            end 
                                                            
            
        end % nofn loop
        
    end % components loop
    
    if count ~= 0 % saving structure to .mat, if any 
        save( strcat('../mat/DRT_',num2str( val ),'_MetricsData','.mat'),'metrics');
        disp('----> metrics .mat file saved.')

        save( strcat('../mat/DRT_',num2str( val ),'_LinRegrData','.mat'),'linregr'); 
        disp('----> regression .mat file saved.')
    else
        disp('----> No components found.');
    end
    
    %clear metrics linregr
    
end % DRT loop

%% ENDINGS
d.printings(d.progStat{2});
dm.deactivateLog;