%% mainDRTVOIGraphData
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 18th, 2015        
%             
%   description: sets up graph and components for the VOI's
%                network based on the DRT field value and stores
%                data structures.

clear all; close all; clc; 
set(0,'DefaultAxesFontSize',18);          

%% 
% load DRT matrix
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;

% directory
matFiles = dir('../mat/DRT_VOI*.mat'); 
numfiles = length(matFiles);

% statistics
drtval = zeros(numfiles,1);
qntComps = zeros(numfiles,1); 
qntGreaterComps = zeros(numfiles,1);

% voxelCoords
bigNetworkCoords = cell(numfiles,1);

% DRT to choose 
drtValue = 9;

% well location
wI = 45; wJ = 68;

% DRT values over the well
drtWell = DRT( wI, wJ,:); 
drtWell = reshape(drtWell, [size(DRT,3) 1]);
Z = find(drtWell == drtValue);

% well coordinate; component id
associatedComps = cell( length(Z),2);
lateralCompVoxels = cell(length(Z),2);

for k = 1:numfiles
    st = load( strcat('../mat/',matFiles(k).name) ); 
    val = st.VOISt.value;     
    drtval(k) = val;
            
    ncomp = length(st.VOISt.compVoxelCoords);
    qntComps(k) = ncomp;
    qntGreaterComps(k) = st.VOISt.compNNodes{1};     
    bigNetworkCoords{k} = st.VOISt.compVoxelCoords{1};

    if val == drtValue 
        
        for i = 1:length(Z)
            
            vxBase = [ wI, wJ, Z(i) ]; % voxel at the well

            idCompsToPush = []; % components to associate
            
            % component table for .dat
            itc = zeros(ncomp,1); 
            ntc = 0*itc;
            
            % component loop
            for idcomp = 1:ncomp
                
                compCoords = st.VOISt.compVoxelCoords{idcomp};                        
                for e = 1:size(compCoords,1)
                    
                    id = strmatch( compCoords(e,:), vxBase );       %#ok<*MATCH2>
                    
                    % base voxel is in the component? if so, pushes the
                    % component and goes to next.
                    if ~isempty(id) 
                        idCompsToPush = [ idCompsToPush; idcomp ];
                        continue;
                    end
                end                                                
                % cvi = st.VOISt.compVoxelInds{idcomp};         
                % plotVoxelGraphComp( DRT,cvi,val,idcomp,0.8);                
               
                itc(idcomp) = idcomp;
                ntc(idcomp) = size(st.VOISt.compVoxelCoords{idcomp},1);
            end
            
        associatedComps{i,1} =  Z(i); % z-depth
        associatedComps{i,2} = idCompsToPush; % component id
        
        % lateral components per z + well voxel
        cp = associatedComps{i,2};        
        lateralCompVoxels{i,1} = Z(i);                  
        cpCoords = st.VOISt.compVoxelCoords{cp};
        izcp = find( cpCoords(:,3) == Z(i) );                        
        lateralCompVoxels{i,2} = cpCoords(izcp,:);
        
        C1 = st.VOISt.compVoxelCoords{1};
        C2 = st.VOISt.compVoxelCoords{2};
        C3 = st.VOISt.compVoxelCoords{3};
                   
        end
        
        %--- save .dat                
        A = [itc,ntc];
                    
        hdr = ['idComp ','nvxComp'];
        fname = strcat('idCompTable_DRT_',num2str(val));            
        dlmwrite(strcat('../dat/subdomain/',fname,'.dat'),hdr,'delimiter','');
        dlmwrite(strcat('../dat/subdomain/',fname,'.dat'),A,'-append','delimiter','\t');            
                        
        % --- plots histogram of components
        % cplot = [1,2,180,369,379]; DRT 14 % components to plot        
        cplot = 1:3; % components to plot        
        lim = length(cplot);
        figure
        hname = strcat('CompHistogram_DRT_',num2str(val));                    
        [icc,ncc,mncc,isol,nisol] = plotCompHist(st.VOISt.value,st.VOISt.compVoxelCoords,lim,0.5,[0.5,100]);
        print(gcf,'-dpdf',strcat('../figs/graphPath/',hname));        
        
        % --- plots the component 
        for q = 1:lim;            
            cvi = st.VOISt.compVoxelInds{ cplot(q) };         
            plotVoxelGraphComp( DRT,cvi,val,cplot(q),0.6,bone,'p', [-37.5,30],'eps');                        
        end
        
    end   
           
    clear VOISt;
end

% gross VOI statistics
M = [ drtval, qntComps, qntGreaterComps ];

% Max col
[ilims,jlims,klims,leni,lenj,lenk,Ws,nvcol,cvcol] = getClusterMaxCol(C2,DRT,drtValue);
