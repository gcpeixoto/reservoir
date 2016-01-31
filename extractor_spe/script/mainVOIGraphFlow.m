%% mainVOIGraphFlow

clear all; close all; clc; 
splshScreenVOIGraphFlow();

%% Load properties 

% load DRT matrix
DRT = replaceInfDRT('../mat/DRT_Field.mat');

% directory
[I,J,K] = setGridBounds(60,220,85); % default

% well
%ic = input('----> Enter I well coordinate: \n');
%jc = input('----> Enter J well coordinate: \n');
ic = 45; jc = 68; 

% DRT 
%drtVal = input('----> Enter DRT: \n');
drtVal = 7;

wellfile = strcat( 'Well_I',num2str(ic),'_J',num2str(jc) );
dbase = strcat( '../mat/',wellfile,'/' );

% metrics data structure
aux = load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_MetricsData.mat') );
metrics = aux.metrics;

aux = load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_LinRegrData.mat') );
linregr = aux.linregr;

% DRT data structure
aux = load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_',wellfile,'.mat') );
VOISt = aux.VOISt;
val = VOISt.value;

% Rerun to avoid assembly
%rerun = input('\n ----> Load .txt files for pressure? [0] no; [1] yes \n');
rerun = 0;
if rerun == 1
    rerun = true;
    dbase = strcat( '../txt/',wellfile,'/' );
    pFiles = dir( strcat(dbase,'/press*.txt') ); 
    numfiles = length(pFiles);
elseif rerun == 0
    rerun = false;
    pFiles = dir( strcat(dbase,'press*.mat') ); 
    numfiles = length(pFiles);
else 
    error('Please, enter 0 or 1.');
end


% loop over pressure files
%for k = 1:numfiles
for k = numfiles
    
    [~,fname,~] = fileparts( strcat(dbase,pFiles(k).name) );        
    flows.well = [ic,jc];
    flows.DRT = drtVal;
    flows.year{k} = str2double( fname(end-3:end) );
    
    % REMARK: remove the 2 first lines out from the CMG .txt file
    
    % load files    
    if rerun    % .txt
        press = load( strcat(dbase,pFiles(k).name),'-ascii' );     
        P = assemble3DPressure( press,I,J,K );    
        %plot3DField(I,J,K,P,'Pressure Field');
        
    else        % .mat
        press = load( strcat(dbase,pFiles(k).name) );     
        P = press.P;
        
    end
    
    % clusters loop    
    %for n = 1:length(metrics.idComp)
    for n = 1
        fprintf('----> Cluster = %d \n',n)

        cvc = VOISt.compVoxelCoords{n};
        cvi = VOISt.compVoxelInds{n};
        %plotPField(cvc,P(cvi),1,drtVal);
                        
        % adjacency matrix
        MadjFlowP = sparse( length(cvc),length(cvc) );
        
        % cluster's elements loop
        for nv = 1:length(cvc)                                      
            [ vc6n, vi6n, p6n ] = getVoxel6Neigh( cvc(nv,:), cvc, cvi, P );                                
                                    
            pc = P( cvi(nv) );  % pressure at the voxel 
            
            %fprintf('-----> Voxel = %d; X = [%d %d %d]; p = %f \n',...
            %         nv,cvc(nv,1),cvc(nv,2),cvc(nv,3),pc);
            
            % finds voxel's 6-neighbours
            for viz = 1:length(vi6n)
                
                    %fprintf('------> Voxel = %d; X = [%d %d %d]; p = %f \n',...
                    %    vi6n(viz),vc6n(viz,1),vc6n(viz,2),vc6n(viz,3),p6n(viz));

                % checking pressure gradient to fill the 
                % adjacency matrix for the directed
                % graph of the flow network
                %
                % 6-neighbours 
                % 
                %                 ____                 
                %                |    |                 
                %                | p1 |                
                %            ____|____|____           
                %           |    |    |    |       
                %           | p2 | pc | p4 |       
                %           |____|____|____|            
                %                |    |                
                %                | p3 |                
                %                |____|
                %
                %
                % Assumes dp = pi - pc i = 1,...,6
                %
                % dp > 0 => flow from pi (source) to pc (target) 
                %        => Madj(i,c) = 1
                %
                % dp < 0 => flow from pc (source) to pi (target) 
                %        => Madj(c,i) = 1 
                % 
                % vi6n(viz) ~= 0 bypasses indices 0
                
                dp = p6n(viz) - pc;                             
                
                if dp > 0                     
                    %disp('-------> dp > 0'); 
                    MadjFlowP( vi6n(viz), nv ) = 1;                        
                    
                elseif dp < 0 
                    %disp('-------> dp < 0'); 
                    MadjFlowP( nv, vi6n(viz) ) = 1;                       
                    
                end                          
                
            end % close 6-neighbour loop
                       

        end % close cluster's element loop  
                                        
        % ----- creating NETWORKX interface                
        
        % max closeness voxel
        closeness{n} = metrics.closenessCentrality{n};
        maxclo{n} = max( closeness{n} );
        imaxclo{n} = find( closeness{n} == maxclo{n} );    
        cvcmaxclo{n} = cvc(imaxclo{n},:);
        
        % min closeness voxel
        minclo{n} = min( closeness{n} );
        iminclo{n} = find( closeness{n} == minclo{n} );
        cvcminclo{n} = cvc(iminclo{n},:);
                        
        % export parameters to networkx
        [capfile,ssfile] = setGraphParams(MadjFlowP, P(cvi),imaxclo{n},iminclo{n} );         
        %maxflow = getFlowData(capfile,ssfile); 
        
        

        %flows.year{k}.cluster{n}.maxflow = maxflow;
    
    end % close cluster's loop
    
    % save .vtk         
    %saveVtkCellCentered( P, fname, 'pressure_cell', 'pressure_point');    

    % save .mat
    if ~rerun
        save( strcat('../mat/',fname,'.mat'), 'P' );                
    end
        
    clear P;

end % end loop of pressure files
