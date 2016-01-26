%% mainVOIGraphFlow

clear all; close all; clc; 

%% Load properties 


aux = load('../mat/PHI.mat');
PHI = aux.PHI;
aux = load('../mat/KN_Field.mat');
KN = aux.KN;

% directory
[I,J,K] = setGridBounds(60,220,85); % default

base_dir = '../txt/press_45_68/';
pFiles = dir( strcat(base_dir,'*.txt') ); 
numfiles = length(pFiles);

drtVal = '7'; % DRT to get

% metrics data structure
aux = load( strcat('../mat/VOI_DRT_',drtVal,'_MetricsData_.mat') );
metrics = aux.metrics; % components with more than 10 voxels

% DRT data structure
aux = load( strcat('../mat/DRT_VOI_',drtVal,'_Well_I45_J68.mat') );
VOISt = aux.VOISt; % all components


% loop over pressure files
for k = 1:numfiles
    
    % REMARK: remove the 2 first lines out from the CMG .txt file
    
    % load files
    press = load( strcat(base_dir,pFiles(k).name),'-ascii' );     
    P = assemble3DPressure( press,I,J,K );    
    %plot3DField(I,J,K,P,'Pressure Field');
    
    % clusters loop    
    %for n = 1:length(metrics.idComp)
    for n = 1

        cvc = VOISt.compVoxelCoords{n};
        cvi = VOISt.compVoxelInds{n};

        MadjFlowP = sparse( length(cvc),length(cvc) );

        % cluster's elements loop
        for nv = 1:length(cvc)
            [ vc6n, vi6n, p6n ] = getVoxel6Neigh( cvc(nv,:), cvc, cvi, P );        

            pc = P( cvi(nv) );  % pressure at the voxel 
            
            % finds voxel' neighbours
            for viz = 1:6

                % checking pressure gradient to fill the 
                % adjacency matrix for the undirected
                % graph of the flow network
                
                if p6n(viz) > pc && p6n(viz) ~= 0 
                    MadjFlowP( vi6n(viz), nv ) = 1;
                elseif p6n(viz) < pc && p6n(viz) ~= 0 
                    MadjFlowP( nv, vi6n(viz) ) = 1;
                end
                                
                
            end % close 6-neighbour loop

        end % close cluster's element loop  
        
        % ----- creating SNAP interface                
        capfile = saveCapacityTable(MadjFlowP, P(cvi) );  
        %! ./../cpp/graphFlow
        %[nodeID,deg,clns,betw] = getMetricsData(capfile);                
        
        % print table of flow network
        

        % save .vtk 
        fname = strcat( 'pressure_time',num2str(k-1) );    
        %saveVtkCellCentered( P, fname, 'pressure_cell', 'pressure_point');    

        % save .mat
        %save( strcat('../mat/',fname,'.mat') );

    end


end % end loop of pressure files
