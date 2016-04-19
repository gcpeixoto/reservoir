%% mainFieldDRTConnections - find DRT connected components for the whole field
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 21st, 2015        
%             
%   description: sets up graph and components for the whole field's
%                network based on the DRT field value and stores
%                data structures.
%
%   requirements:
%        - pre-computed .mat files
%        - Gibbon toolbox functions
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

% voxel connectivity
% REMARK: 26-neigh is invalid for CMG (no flow; finite volume approach)
dv = setNeighDist('6'); 

csv_flag = false; % to save .csv for whole field 
vtk_flag = false; % to export .vtk of DRT field 

% minimum number of voxels desired per component to save .csv table
nofsv = 1000; 

%% LOAD .MAT FILES

[PHI,KX,KY,KZ,KN,PHIZ,RQI,FZI,DRT] = loadMatFiles;

%% REAPING DRTs  
drt = sort( unique( DRT(:) ) );
 
%{
    tabulate (Statistics data)
    DRT value | frequency | percentage
%}
tab = tabulate( DRT(:) );
countDRT = tab(:,2);

% ignores DRTs whose frequency is below 100
Ig = find( countDRT(:) < 100 );
drtIgnored = drt(Ig);
C = find( countDRT(:) > 100 );
drt = drt(C); 
drt = drt(drt>0); % removes 0

% file header used in the loop
head = {'i,'; 'j,'; 'k,'; 'phi_e,'; ...
        'kx,'; 'ky,'; 'kz,'; 'kn,'; ...
        'phiz,'; 'RQI,'; 'FZI,'; 'LogPhiz,'; 'LogRQI'}'; % transposed!
txt = sprintf('%s\t',head{:});
txt(end)='';

%%% --------------------- Sweeping DRTs

for m = 1:length(drt)
               
    fprintf('----> Sweeping field: DRT = %d... \n',drt(m));
    
    % filtering grid to capture voxels with a specific DRT
    [ VBDRT, coordsDRT, indz ] = maskVoxels( DRT,drt(m) ); 
                
    %{
        Adjacency matrix
        ================
                    
        Strategy:        
                    
            - Find the entries whose distance is less than sqrt(3), which
              stands for the 6-voxel (or 26-voxel) neighbourhood
    
            - Set up sparse adjacency matrix from the entries found to be
              able to create an undirected graph
    
            - Store graph edge list matrix
              
    %}
            
    indIJ = [];
    for i = 1:size(coordsDRT,1)
        %fprintf('----> i = %d... \n',size(coordsDRT,1)-i); 
        for j = i:size(coordsDRT,1)
                                                            
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
    aux = [ indIJ(:,2) indIJ(:,1) ]; % reverse edges [ j i ]
    indIJ = [ indIJ; aux ]; % filling
    
    if isempty(indIJ)
        fprintf('----> No connections found for DRT = %d... \n',drt(m));
        continue;
    end
    
    disp('----> Computing adjacency matrix...');            
    % creates adjacency matrix n x n by marking 1 for connected nodes
    MDadj = sparse( indIJ(:,1),indIJ(:,2),1,size(coordsDRT,1),size(coordsDRT,1) ); 
            
    %{ 
        Finding network components 
        ==========================     
    
        - Uses the function 'networkComponents', by Daniel Larremore 
        @MathWorks central
    
        - From adjacency matrix, computes all the components 
          (isolated + connected) for the specific DRT network
    
        - The function is based on the connected component algorithm
          for graphs (e.g., see 
         http://people.sc.fsu.edu/~jburkardt/classes/asa_2011/asa_2011_graphs.pdf)
    
       
        DRT Data structure
        ==================
    
                drtSt
                  |           
                  |- value  
                  |
                  |- all{ ... }                 }
                  .                             }
                  .                             } Global family ( all voxels with such DRT )
                  .                             }
                  |                             } [ SEVERAL DATA ]
                  |
                  |- comp{ ... }                }        
                  |    |                        }
                  |    |- comp{ ... }{idcomp1}  }
                  |    |- comp{ ... }{idcomp2}  } 
                  .    .                        } Cluster Family ( K components of variable voxels per DRT )
                  .    .                        }
                  .    .                        }
                  |    |                        }
                  |    |- comp{ ... }{idcompK}  } [ SEVERAL DATA ]
                                                                  
    %}      
    disp('----> Finding network components...');
    [ncomp,compSizes,members] = networkComponents(MDadj);
            
    % global  
    drtSt.value = drt(m);                     % DRT value
    drtSt.allAdjMatrix = MDadj;               % graph adjacency matrix
    drtSt.allAdjEdgeList = edgeList;          % graph edge list
    drtSt.allVoxelCoords = coordsDRT;         % voxel coordinates (i,j,k)     
    drtSt.allVoxelInds = indz;                % voxel linear indices
    drtSt.allPHI = PHI( indz );               % effective porosity
    drtSt.allKX = KX( indz );                 % permeability-x
    drtSt.allKY = KY( indz );                 % permeability-y
    drtSt.allKZ = KZ( indz );                 % permeability-z
    drtSt.allKN = KN( indz );                 % permeability-norm
    drtSt.allPHIZ = PHIZ( indz );             % normalised porosity
    drtSt.allRQI = RQI( indz );               % RQI
    drtSt.allFZI = FZI( indz );               % FZI
    drtSt.allLogPHIZ = log10( PHIZ( indz ) ); % as is
    drtSt.allLogRQI = log10( RQI( indz ) );   % as is
    drtSt.allNComps = ncomp;                  % number of components in the graph

    % mount matrix to export
    mat = [ coordsDRT(:,1) ... 
            coordsDRT(:,2) ...
            coordsDRT(:,3) ...
            PHI( indz )    ...
            KX( indz )     ...
            KY( indz )     ...
            KZ( indz )     ...
            KN( indz )     ...
            PHIZ( indz )   ...
            RQI( indz )    ...
            FZI( indz )    ...
     log10( PHIZ( indz ) ) ...
     log10(  RQI( indz ) ) ];
            
    if csv_flag == true
        
        % preparing csv file
        fname = strcat('../csv/Table-Field','_DRT_',num2str( drt(m) ),'.csv');    
        dlmwrite(fname,txt,'');        
    
        % append matrix
        dlmwrite(fname,mat,'-append');
        fprintf('----> Table-Field_DRT_%s.csv file saved. \n',...
                num2str(drt(m)));
        disp('----> .csv file saved - global.');
    end
    
    % cluster (each component)
    for idcomp = 1:ncomp
        
        fprintf('--------> Computing component = %d... \n',idcomp);
        
        aux = coordsDRT( members{idcomp},: );
        drtSt.compVoxelCoords{idcomp} = aux; 
        drtSt.compVoxelInds{idcomp} = indz( members{idcomp} );        
        drtSt.compNNodes{idcomp} = compSizes(idcomp);        
        drtSt.compPHI{idcomp} = PHI( indz( members{idcomp} ) );
        drtSt.compKX{idcomp} = KX( indz( members{idcomp} ) );
        drtSt.compKY{idcomp} = KY( indz( members{idcomp} ) );
        drtSt.compKZ{idcomp} = KZ( indz( members{idcomp} ) );
        drtSt.compKN{idcomp} = KN( indz( members{idcomp} ) );
        drtSt.compPHIZ{idcomp} = PHIZ( indz( members{idcomp} ) );
        drtSt.compRQI{idcomp} = RQI( indz( members{idcomp} ) );
        drtSt.compFZI{idcomp} = FZI( indz( members{idcomp} ) );
        drtSt.compLogPHIZ{idcomp} = log10( PHIZ( indz( members{idcomp} ) ) );
        drtSt.compLogRQI{idcomp} = log10( RQI( indz( members{idcomp} ) ) ); 
                  
        % matrix to export
        mat = [ aux(:,1)                         ... 
                aux(:,2)                         ...
                aux(:,3)                         ...
                PHI( indz( members{idcomp} ) )   ...
                 KX( indz( members{idcomp} ) )   ...
                 KY( indz( members{idcomp} ) )   ...
                 KZ( indz( members{idcomp} ) )   ...
                 KN( indz( members{idcomp} ) )   ...
               PHIZ( indz( members{idcomp} ) )   ...
                RQI( indz( members{idcomp} ) )   ...
                FZI( indz( members{idcomp} ) )   ...
        log10( PHIZ( indz( members{idcomp} ) ) ) ...
        log10(  RQI( indz( members{idcomp} ) ) ) ];
        
        if (csv_flag == true) && (compSizes(idcomp) >= nofsv)              
            % preparing csv file
            fname = strcat('../csv/Table-Cluster_',num2str(idcomp),'_DRT_',num2str( drt(m) ),'.csv');        
            dlmwrite(fname,txt,'');       
            
            % append matrix 
            dlmwrite(fname,mat,'-append'); 
            fprintf('----> Table-Cluster_%s_DRT_%s.csv file saved. \n',...
                num2str(idcomp),num2str(drt(m)));            
        end                        
    end              

    % saving structure to .mat     
    save( strcat('../mat/DRT_',num2str( drt(m) ),'.mat'),'drtSt'); % saving
    fprintf('----> DRT_%s.mat file saved. \n',num2str(drt(m)));
    
    clear drtSt;    % frees to recompute
    
end

% export to VTK
if vtk_flag == true
    saveVtkCellCentered(DRT,'DRT_Field','DRT_Voxel','DRT_Point');
end

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;