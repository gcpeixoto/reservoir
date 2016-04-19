%% mainVOIDRTGraphMetrics
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 17th, 2015        
%             
%   description: compute metrics and regression analysis for the VOI
%                networks based on the DRT field value.
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
d.VOIgraphDataDependency;

%% LOAD FILES

[~,~,~,~,~,PHIZ,RQI,FZI,DRT] = loadMatFiles;

% well 
ic = 45; jc = 68;

dbase = strcat( '../mat/Well_I',num2str(ic),'_J',num2str(jc),'/' );
matFiles = dir( strcat(dbase,'VOI_DRT*.mat') ); 
numfiles = length(matFiles);
    
% sweeping DRTs
for k = 1:numfiles 
    
    st = load( strcat(dbase,matFiles(k).name) ); 
    val = st.VOISt.value; 
    
    fprintf('----> Sweeping DRT: %d... \n',val);
    
    avc = st.VOISt.allVoxelCoords;
    ncomps = st.VOISt.allNComps;    
    Madj = st.VOISt.allAdjMatrix;    
    
    metrics.drtValue = val;
    linregr.drtValue = val;
    
    count = 0;
    for idComp = 1:ncomps        
        
        cnn = st.VOISt.compNNodes{idComp};        
                        
        if cnn > 10 % components with more than 10 nodes
            
            cvc = st.VOISt.compVoxelCoords{idComp};
            cvi = st.VOISt.compVoxelInds{idComp};            
  
            % performs linear regression
            logPHIZ = log10( PHIZ(cvi) );
            logRQI  = log10( RQI(cvi) );
            [ R, m, b ] = regression( logPHIZ, logRQI, 'one' );
                        
            
            % linear regression criteria            
%            if ( m >= 0.95 && m <= 1.05 ) && ( R*R >= 0.9 && R*R <= 1.0 )
                          
                fprintf('----> Good component found: %d. Computing subgraph... \n',idComp);
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
                                
                % ------------- @MIT Strat. Eng. VERY SLOW! 
                %disp('----> Computing degree centrality...');
                %[deg,~,~] = degrees(MadjComp);            % degree centrality                
                
                %disp('----> Computing closeness centrality...');
                %clns = closeness(MadjComp);               % closeness centrality   
                
                %disp('----> Computing betweeness centrality...');
                %betw = node_betweenness_slow(MadjComp);   % betweeness centrality                
                %betw = betw';
                
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
                metrics.centerVoxelCoords{count} = ivC;
                metrics.adjMatrix{count} = MadjComp;
                                
                linregr.idComp{count} = idComp;
                linregr.Pearson{count} = R*R;
                linregr.slope{count} = m;
                linregr.offset{count} = b;
                linregr.logPHIZ{count} = logPHIZ;
                linregr.logRQI{count} = logRQI;                                
                
%            end % regression loop
            
        end % components with > 10 loop
        
    end % components loop
    
    if count ~= 0 % saving structure to .mat, if any 
        save( strcat(dbase,'VOI_DRT_',num2str( val ),'_MetricsData','.mat'),'metrics');
        disp('----> metrics .mat file saved.')

        save( strcat(dbase,'VOI_DRT_',num2str( val ),'_LinRegrData','.mat'),'linregr'); 
        disp('----> regression .mat file saved.')
    else
        disp('----> No components found.');
    end
    
    clear VOISt;
    
end % DRT loop

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;