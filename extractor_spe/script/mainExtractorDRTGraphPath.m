%% extractorDRTGraphPath
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Oct 21st, 2015        
%             
%   description: sets up graphs for the reservoir network based on 
%                the DRT field value.
%
%   requirements:
%        - pre-computed .mat files
%        - Gibbon toolbox functions
%        - Matlab third-party additional functions
%

%% DEFAULTS
clear all; close all; clc; format long;
delete('../log/extractorDRTGraphPath.log');
diary('../log/extractorDRTGraphPath.log');
diary on
%profile on 

setOptions;

%% LOAD FILES

% file paths
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

[PHI,KX,KY,KZ] = loadMatFiles(phiname,kxname,kyname,kzname);

%% COMPUTATION OF PROPERTIES

KN = sqrt( KX.^2 + KY.^2 + KZ.^2 );
PHIZ = PHI./(1.0 - PHI);
RQI = 0.0314*sqrt( KN./PHI );
FZI = RQI./PHIZ;

[ VB, coords, inds ] = maskVoxels( FZI,Inf );
if ~isempty( coords ) 
    for e = 1:size(coords,1)
        FZI( coords(e,1), coords(e,2), coords(e,3) ) = 0.0;
    end
end

DRT = round( 2*log( FZI ) + 10.6 );
save('../mat/DRT_Field.mat','DRT'); % saving

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
drt = drt(2:end); % removes -Inf (due to phi = 0 values)


% file header used in the loop
head = {'i,'; 'j,'; 'k,'; 'phi_e,'; 'kx,'; 'ky,'; 'kz,'; 'kn,'; 'phiz,'; 'RQI,'; 'FZI,'; 'LogPhiz,'; 'LogRQI'}; 
head = head';
txt=sprintf('%s\t',head{:});
txt(end)='';

%%% --------------------- Sweeping DRTs

for m = 1:length(drt)
               
    fprintf('----> Computing DRT = %d... \n',drt(m));
    
    % filtering grid to capture voxels with a specific DRT
    [ VBDRT, coordsDRT, indz ] = maskVoxels( DRT,drt(m) ); 
                
    %{
        Adjacency matrix
        ================
                    
        Strategy:        
                    
            - Find the entries whose distance is less than sqrt(3), which
              stands for the 26-voxel neighbourhood
    
            - Set up sparse adjacency matrix from the entries found to be
              able to create an undirected graph
    
            - Store graph edge list matrix
              
    %}
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
        
    disp('----> Setting up structure - global ...');
    % global  
    drtSt.value = drt(m);                    % DRT value
    drtSt.allAdjMatrix = MDadj;              % graph adjacency matrix
    drtSt.allAdjEdgeList = edgeList;         % graph edge list
    drtSt.allVoxelCoords = coordsDRT;        % voxel coordinates (i,j,k)     
    drtSt.allVoxelInds = indz;               % voxel linear indices
    drtSt.allPHI = PHI( indz );              % effective porosity
    drtSt.allKX = KX( indz );                % permeability-x
    drtSt.allKY = KY( indz );                % permeability-y
    drtSt.allKZ = KZ( indz );                % permeability-z
    drtSt.allKN = KN( indz );                % permeability-norm
    drtSt.allPHIZ = PHIZ( indz );            % normalised porosity
    drtSt.allRQI = RQI( indz );              % RQI
    drtSt.allFZI = FZI( indz );              % FZI
    drtSt.allLogPHIZ = log( PHIZ( indz ) );  % as is
    drtSt.allLogRQI = log( RQI( indz ) );    % as is
    drtSt.allNComps = ncomp;                 % number of components in the graph

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
       log( PHIZ( indz ) ) ...
       log(  RQI( indz ) ) ];
        
    % preparing csv file
    fname1 = strcat('../csv/graphpath/GraphDataAll','_DRT_',num2str( drt(m) ),'.csv');    
    dlmwrite(fname1,txt,'');        
    
    % append matrix
    dlmwrite(fname1,mat,'-append'); 
    disp('----> csv file saved - global.');
          
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
        drtSt.compLogPHIZ{idcomp} = log( PHIZ( indz( members{idcomp} ) ) );
        drtSt.compLogRQI{idcomp} = log( RQI( indz( members{idcomp} ) ) ); 
                  
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
          log( PHIZ( indz( members{idcomp} ) ) ) ...
          log(  RQI( indz( members{idcomp} ) ) ) ];
                      
        % getting only the most connected component (first one in the list)
        if idcomp == 1 
            disp('----> csv file saved - component.');
            % preparing csv file
            fname2 = strcat('../csv/graphpath/GraphDataComp_',num2str(idcomp),'_DRT_',num2str( drt(m) ),'.csv');        
            dlmwrite(fname2,txt,'');        
         
            % append matrix 
            dlmwrite(fname2,mat,'-append'); 
        end
        
                
    end              

    %plotVoxelGraphComp( DRT,drtSt.compVoxelInds{1},drtSt.value,1,0.8);        

    % saving structure to .mat     
    save( strcat('../mat/DRT_',num2str( drt(m) ),'.mat'),'drtSt'); % saving
    disp('----> .mat file saved.')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING
%     pltfield = false; % plot DRT = drt for all the field?
%     % only graph and paths with more than 'nvpath' voxels
%     nvpath = 2;
%     if size(path,1) > nvpath
%         
%         %------------- GRAPH PATH
%         figure     
%         path1 = path(:,1);
%         path2 = path(:,2);
%         path3 = path(:,3);        
%         plot3(path2,path1,path3,'LineWidth',1.5, 'Color',[0 .4 0.54], ...
%             'Marker','o', 'MarkerSize',2.5, ...
%             'MarkerFaceColor','k', 'MarkerEdgeColor','k')
% 
%         axis equal, axis vis3d, grid on, view(3)
%         title( strcat( 'Path - DRT ', num2str( drt(m) ) ) );
%         xlabel('J');
%         ylabel('I');
%         zlabel('K');
%         pt1 = sort( unique( path1 ) );
%         pt2 = sort( unique( path2 ) );
%         pt3 = sort( unique( path3 ) );
%         pt1 = num2str(pt1); pt2 = num2str(pt2); pt3 = num2str(pt3);
%         pt1 = cellstr(pt1); pt2 = cellstr(pt2); pt3 = cellstr(pt3);
%         
%         set( gca,'XTick',min(path2):1:max(path2),'XTickLabel',pt2, ...
%                  'YTick',min(path1):1:max(path1),'YTickLabel',pt1, ... 
%                  'ZTick',min(path3):1:max(path3),'ZTickLabel',pt3 );
%         
%         print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
%                 strcat('Path_DRT_',num2str( drt(m) ) ) ) );   
%         
%             
%         %------------- VOXEL PATH
%         figure 
%         [ ~, ~, indz ] = maskVoxels( DRT,drt(m) );
%         [F,V,C]=ind2patch(idpath,DRT,'v');
%         
%         hold on;
%         patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.8);
%         axis equal; view(3); axis tight; axis vis3d; grid off;        
%         title( strcat( 'Path - DRT ', num2str( drt(m) ) ) );
%         xlabel('J');
%         ylabel('I');
%         zlabel('K');        
%         
%         print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
%                 strcat('Voxel_Path_DRT_',num2str( drt(m) ) ) ) ); 
%             
%         if pltfield == true
%             plotDRTField(DRT,val);
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%profsave(profile('info'),'../log/profile_report')
close all
diary off