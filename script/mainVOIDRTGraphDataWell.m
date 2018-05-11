%% mainVOIDRTGraphDataWell - compute statistics of volumetrics
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 18th, 2015        
%             
%   description: computes statistics of volumetrics for clusters 
%                that intersect the well of interest (central well
%                of the VOI) inside the chosen VOI.
%                In the paper, these are the well-intersection clusters.
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
d.VOIConnectionsDependency;

%% INPUTS

ic = 45; jc = 68;   % center well coordinates
drtValue = 13;      % DRT to choose 

%% LOAD FILES

[~,~,~,~,~,~,~,~,DRT] = loadMatFiles;

% directory
dbase = strcat( '../mat/Well_I',num2str(ic),'_J',num2str(jc),'/' );

matFiles = dir( strcat(dbase,'VOI_DRT*.mat') ); 
matFiles = checkMetricsFiles(matFiles,dbase); % required because 'VOISt'
numfiles = length(matFiles);

% DRT values over the well
drtWell = DRT(ic,jc,:); drtWell = reshape(drtWell, [size(DRT,3), 1]);

% DRT depths 
Z = find(drtWell == drtValue); 

% well coordinate; component id
associatedComps = cell( length(Z),2);
lateralCompVoxels = cell(length(Z),2);

% dir
cld = '../csv/clusterDataWell/';
if exist(cld,'dir') ~= 7; mkdir(cld); end   

% sweep DRT clusters
for k = 1:numfiles
    
    load( strcat(dbase,matFiles(k).name), 'VOISt' ); 
    val = VOISt.value;               
    ncomp = length(VOISt.compVoxelCoords);    
    
    if val == drtValue 
        
        for i = 1:length(Z)         % sweep by depths
            
            vxBase = [ic,jc,Z(i)];  % voxel at the well
            idCompsToPush = [];     % components to associate
            
            % component table for file
            itc = zeros(ncomp,1);   % cluster index
            ntc = 0*itc;            % number of voxels                                     
                        
            for idcomp = 1:ncomp    % component loop                                

                compCoords = VOISt.compVoxelCoords{idcomp};                                        
                
                for e = 1:size(compCoords,1)                                        
                
                    id = strmatch(compCoords(e,:),vxBase);  %#ok<*MATCH2>                    
                    % base voxel is in the component? if so, pushes the
                    % component and goes to next.
                    if ~isempty(id) 
                        idCompsToPush = [idCompsToPush;idcomp];
                        continue;
                    end
                end   
                                                               
                itc(idcomp) = idcomp;
                ntc(idcomp) = size(VOISt.compVoxelCoords{idcomp},1);                
                
            end
                                    
        %{            
                well
                 |
                 _
                |_|
           _ _ _|_|_ _
          |d|d|d|d|d|d| <--- lateral components vi for z=zi with DRT = d
                |_|          (includes well-intersection voxel)
                |_|            
                 .
                 .
                 .
        
            - associatedComps: 
                                    
            [z=z1] [vci] (for z=z1) ---- cluster index vci 
            [z=z2] [vcj] (for z=z2) ---- cluster index vcj 
              .            .             .            
              .            .             .            
              .            .             .            
            [z=zk] [vck] (for z=zk) ---- cluster index vck 
            
            
            - lateralCompVoxels:
        
            [z=z1] [v1,v2,v3,...vn1] (for z=z1) ---- vn1 voxels
            [z=z2] [v1,v2,v3,...vn2] (for z=z2) ---- vn2 voxels
              .            .             .            .
              .            .             .            .
              .            .             .            .
            [z=zk] [v1,v2,v3,...vnk] (for z=zk) ---- vnk voxels
        
        %}
            
        associatedComps{i,1} =  Z(i);           
        associatedComps{i,2} = idCompsToPush;   
        
        cp = associatedComps{i,2};        
        lateralCompVoxels{i,1} = Z(i);                  
        cpCoords = VOISt.compVoxelCoords{cp};
        izcp = find( cpCoords(:,3) == Z(i) );                        
        lateralCompVoxels{i,2} = cpCoords(izcp,:);
        
        %C1 = VOISt.compVoxelCoords{1};        
                   
        end
        
        % mounting table to file 
        nvres = VOISt.allNVoxels; % number of voxels in the reservoir
        nf = cumsum(ntc);              
        vcf = ntc/nf(end);        % volume fraction: cluster-to-DRT family
        vcr = ntc/nvres;          % volume fraction: cluster-to-reservoir
        TBL = [itc,ntc,vcf,vcr];  % volumetrics table
        
        % write to .csv
        hdr = {'idComp,';'nComp,';'nCompi/nComp,';'nCompi/nRes'};
        hdr = hdr';            
        fname = strcat('ClusterTable_DRT_',num2str(val));                                   
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(cld,fname,'.csv'),TBL,'-append','delimiter',',');  
        
        % mounting table to file (well-intersection clusters)
        iac = zeros(size(associatedComps,1),1);
        eac = 0*iac;
        fac = 0*iac;
        rac = 0*iac;
        wzc = 0*iac;
        for nac = 1:size(associatedComps,1)
            iac(nac) = associatedComps{nac,2};
            eac(nac) = ntc(associatedComps{nac,2});
            fac(nac) = vcf(associatedComps{nac,2});
            rac(nac) = vcr(associatedComps{nac,2});
            wzc(nac) = associatedComps{nac,1};
        end           
        TBL = [iac,eac,fac,rac,wzc]; 
            
        % write to .csv
        hdr = {'idComp,';'nComp,';'nCompi/nComp,';'nCompi/nRes,';'intersec_z'};
        hdr = hdr';            
        fname = strcat('IntersectionCompTable_DRT_',num2str(val));            
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(cld,fname,'.csv'),TBL,'-append','delimiter',',');       
                        
        % --- plots histogram of components                
        cplot = 1:4; % components to plot        
        lim = length(cplot);
        hname = strcat('CompHistogram_DRT_',num2str(val));                    
        %figure        
        %[icc,ncc,mncc,isol,nisol] = plotCompHist(VOISt.value,VOISt.compVoxelCoords,lim,0.5,[0.5,100]);
        %print(gcf,'-dpdf',strcat('../figs/',hname));        
        
        % --- plots the component 
        %for q = 1:lim;            
        %    cvi = VOISt.compVoxelInds{ cplot(q) };         
        %    plotVoxelGraphComp( DRT,cvi,val,cplot(q),0.6,bone,'p', [-37.5,30],'pdf');                        
        %end
        
    end   
           
    clear VOISt;
end

% Max col
%[ilims1,jlims1,klims1,leni1,lenj1,lenk1,Ws1,nvcol1,cvcol1] = getClusterMaxCol(C1,DRT,drtValue);

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
