%% mainDRTGraphMetrics
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 27st, 2015        
%             
%   description: compute metrics and regression analysis for the reservoir 
%                networks based on the DRT field value.
%
%   requirements:
%        - pre-computed .mat files
%        - Matlab third-party additional functions
%

%% DEFAULTS
clear all; close all; clc; format long;
delete('../log/DRTgraphMetrics.log');
diary('../log/DRTgraphMetrics.log');
diary on

disp('---- E X E C U T I N G   G R A P H   M E T R I C S ----');

aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;

matFiles = dir('../mat/DRT_*.mat'); 
numfiles = length(matFiles);

% azim, elev angles for 3D plot (by test)
% az = -37.5*ones(1,numfiles);
% el = 30*ones(1,numfiles);
% DRT 4;
% az(4) = - 142;
% el(4) = 22;

% Marker{Face,Edge}Color (scatter)
%mfcolor = [0.5 0.1 0.3];
%mecolor = [0.1 0.1 1.0];
    
% sweeping DRTs

for k = 1:numfiles 
    
    st = load( strcat('../mat/',matFiles(k).name) ); 
    val = st.drtSt.value;    
    fprintf('----> Sweeping DRT: %d... \n',val);
    
    avc = st.drtSt.allVoxelCoords;
    ncomps = st.drtSt.allNComps;    
    Madj = st.drtSt.allAdjMatrix;    
    
    metrics.drtValue = val;
    linregr.drtValue = val;
    
    count = 0;
    for idComp = 1:ncomps        
        
        cnn = st.drtSt.compNNodes{idComp};        
                        
        if cnn > 10 % components with more than 10 nodes
            
            cvc = st.drtSt.compVoxelCoords{idComp};
            cvi = st.drtSt.compVoxelInds{idComp};            
  
            % performs linear regression
            logPHIZ = st.drtSt.compLogPHIZ{idComp};
            logRQI  = st.drtSt.compLogRQI{idComp};
            [ R, m, b ] = regression( logPHIZ, logRQI, 'one' );
                        
            
            % linear regression criteria            
            if ( m >= 0.95 && m <= 1.05 ) && ( R*R >= 0.9 && R*R <= 1.0 )
                          
                fprintf('----> Good component found: %d. Computing data... \n',idComp);
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

                % @MIT Strat. Eng. 
                [deg,~,~] = degrees(MadjComp);            % degree centrality
                clns = closeness(MadjComp);               % closeness centrality  
                betw = node_betweenness_slow(MadjComp);   % betweeness centrality
                betw = betw';
                
                maxC = max(clns);                         % max closeness = min farness
                iC = find( clns == maxC );                % network closer nodes
                ivC = avc( v(iC),: );                     % global voxel coordinates
                                
                % store good components         
                metrics.idComp{count} = idComp;                
                metrics.degreeCentrality{count} = deg;
                metrics.closenessCentrality{count} = clns;
                metrics.betweenessCentrality{count} = betw;
                                
                linregr.idComp{count} = idComp;
                linregr.Pearson{count} = R*R;
                linregr.slope{count} = m;
                linregr.offset{count} = b;
                linregr.logPHIZ{count} = logPHIZ;
                linregr.logRQI{count} = logRQI;
                
                %%%%%%%%%%%%%%%%%%%%%%%% plotting to highlight closeness points         
                %     figure      
                %     [F,V,C]=ind2patch(cvi,DRT,'v');    
                %     hold on;
                %     patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.1);
                %     axis equal; view( [ az(val), el(val) ] ); %view(3);     
                %     axis tight; axis vis3d; grid off;            
                %     
                %     title( strcat( 'Voxel Component: ',num2str(idComp),'- DRT ', num2str( val ) ) );
                %     xlabel('J');
                %     ylabel('I');
                %     zlabel('K');
                %     colormap('gray');
                % 
                %     % scatter plot     
                %     hold on        
                %     scatter3( ivC(:,2), ivC(:,1), ivC(:,3), 300, 'fill', 'MarkerFaceColor', mfcolor,...
                %                                                         'MarkerEdgeColor', mecolor ) ;             

                %     print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
                %           strcat('Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ) ) ) );  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end % regression loop
            
        end % components with > 10 loop
        
    end % components loop
    
    if count ~= 0 % saving structure to .mat, if any 
        save( strcat('../mat/DRT_',num2str( val ),'_MetricsData_','.mat'),'metrics');
        disp('----> metrics .mat file saved.')

        save( strcat('../mat/DRT_',num2str( val ),'_LinRegrData_','.mat'),'linregr'); 
        disp('----> regression .mat file saved.')
    else
        disp('----> No components found.');
    end
    
end % DRT loop

%close all
diary off
disp('---- N O R M A L   T E R M I N A T I O N ----');