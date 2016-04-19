%% mainVOIDRTGraphData
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

% classes
dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEDependency;
d.graphDataDependency;

%% 
% load DRT
[~,~,~,~,~,~,~,~,DRT] = loadMatFiles;

% well
ic = 45; jc = 68;
 
% directory
dbase = strcat( '../mat/Well_I',num2str(ic),'_J',num2str(jc),'/' );
matFiles = dir( strcat(dbase,'VOI_DRT_*_Well*','.mat') ); 
numfiles = length(matFiles);

% statistics
drtval = zeros(numfiles,1);
qntComps = zeros(numfiles,1); 
qntGreaterComps = zeros(numfiles,1);

% voxelCoords
bigNetworkCoords = cell(numfiles,1);

% number of voxels at reservoir 
% TODO automatize
nvres = 29*29*85;

% DRT to choose 
drtValue = 14;

% DRT values over the well
drtWell = DRT( ic, jc,:); 
drtWell = reshape(drtWell, [size(DRT,3) 1]);
Z = find(drtWell == drtValue);

% well coordinate; component id
associatedComps = cell( length(Z),2);
lateralCompVoxels = cell(length(Z),2);

for k = 1:numfiles
    st = load( strcat(dbase,matFiles(k).name) ); 
    val = st.VOISt.value;     
    drtval(k) = val;
            
    ncomp = length(st.VOISt.compVoxelCoords);
    qntComps(k) = ncomp;
    qntGreaterComps(k) = st.VOISt.compNNodes{1};     
    bigNetworkCoords{k} = st.VOISt.compVoxelCoords{1};    
    
    if val == drtValue 
        
        for i = 1:length(Z)
            
            vxBase = [ ic, jc, Z(i) ]; % voxel at the well

            idCompsToPush = []; % components to associate
            
            % component table for .dat
            itc = zeros(ncomp,1);   % cluster index
            ntc = 0*itc;            % number of voxels                                     
            
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
        
        %C1 = st.VOISt.compVoxelCoords{1};
        %C2 = st.VOISt.compVoxelCoords{2};
        %C3 = st.VOISt.compVoxelCoords{3};
        %C4 = st.VOISt.compVoxelCoords{4};
                   
        end
        
        %--- save .dat
        nf = cumsum(ntc);              
        vcf = ntc/nf(end);  % volume fraction: cluster-to-DRT family
        vcr = ntc/nvres;    % volume fraction: cluster-to-reservoir
        A = [itc,ntc,vcf,vcr];
        
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
        B = [iac,eac,fac,rac,wzc]; % associated comp table
                    
        hdr = {'idComp,';'nComp,';'nCompi/nComp,';'nCompi/nRes'};
        hdr = hdr';            
        fname = strcat('idCompTable_DRT_',num2str(val));            
        
        cld = '../dat/clusterData/';
        if exist(cld,'dir') ~= 7; mkdir(cld); end   % create            
        dlmwrite(strcat(cld,fname,'.dat'),hdr,'delimiter','');
        dlmwrite(strcat(cld,fname,'.dat'),A,'-append','delimiter','\t');            
        dlmwrite(strcat(cld,fname,'.dat'),B,'-append','roffset',1,'delimiter','\t');            
                        
        % --- plots histogram of components                
        cplot = 1:4; % components to plot        
        lim = length(cplot);
        figure
        hname = strcat('CompHistogram_DRT_',num2str(val));                    
        %[icc,ncc,mncc,isol,nisol] = plotCompHist(st.VOISt.value,st.VOISt.compVoxelCoords,lim,0.5,[0.5,100]);
        %print(gcf,'-dpdf',strcat('../figs/',hname));        
        
        % --- plots the component 
        for q = 1:lim;            
            cvi = st.VOISt.compVoxelInds{ cplot(q) };         
            plotVoxelGraphComp( DRT,cvi,val,cplot(q),0.6,bone,'p', [-37.5,30],'pdf');                        
        end
        
    end   
           
    clear VOISt;
end

% gross VOI statistics
M = [ drtval, qntComps, qntGreaterComps ];

% Max col
%[ilims1,jlims1,klims1,leni1,lenj1,lenk1,Ws1,nvcol1,cvcol1] = getClusterMaxCol(C1,DRT,drtValue);
%[ilims2,jlims2,klims2,leni2,lenj2,lenk2,Ws2,nvcol2,cvcol2] = getClusterMaxCol(C2,DRT,drtValue);
%[ilims3,jlims3,klims3,leni3,lenj3,lenk3,Ws3,nvcol3,cvcol3] = getClusterMaxCol(C3,DRT,drtValue);
%[ilims4,jlims4,klims4,leni4,lenj4,lenk4,Ws4,nvcol4,cvcol4] = getClusterMaxCol(C4,DRT,drtValue);

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
