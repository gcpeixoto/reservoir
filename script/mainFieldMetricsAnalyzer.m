    %% mainFieldMetricsAnalyzer - analyzes metrics in the field
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Aug 29th, 2016        
%             
%   description: finds points of minimum and maximum centralities, plots 
%                the fields and export data to files


  
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
d.fieldGraphMetricsDependency;

%% INPUT  

%{
    How to use:

        - drtVal: array of DRT values to get. Unique or list.        
        e.g.: 
                drtVal = 5; loads only DRT = 5;
                drtVal = [2,5]; loads DRT = 2 and DRT = 5;
                drtVal = 4:10; loads DRTs from 4 to 10.


        - nc: limit number of components to go through per cluster.              
        e.g.:
                nc = []; loads ALL the clusters for a unique DRT;
                nc = 10; loads AT MOST 10 clusters per DRT; 
                nc = 100; loads AT MOST 100 clusters per DRT;          
                      
                If 'nc' is greater than the number of clusters, 
                the code will sweep ALL the clusters.
                        
%}

drtVal = 10:11;
nc = []; 

%% LOAD FILES

% dir  checking
cld = '../csv/fieldData/';
if exist(cld,'dir') ~= 7; mkdir(cld); end   

% load 
[~,~,~,~,KN,~,~,~,DRT] = loadMatFiles;

%% LOOP OVER DRT VALUES

for dv = 1:length(drtVal)
            
    % metrics data structure
    load( strcat('../mat/DRT_',num2str(drtVal(dv)),'_MetricsData.mat'),'metrics' );
    load( strcat('../mat/DRT_',num2str(drtVal(dv)),'_LinRegrData.mat'),'linregr' );

    % DRT data structure
    load( strcat('../mat/DRT_',num2str(drtVal(dv)),'.mat'),'drtSt' );

    ncomp = numel(metrics.idComp); % number of components
    
    if ~isempty(nc) && nc <= ncomp, ncomp = nc; end        
    
    for n = 1:ncomp % cluster loop                   

        compVoxelInds{n} = drtSt.compVoxelInds{ metrics.idComp{n} };        
        compVoxelCoords{n} = drtSt.compVoxelCoords{ metrics.idComp{n} };    
        adjMat{n} = metrics.adjMatrix{n};

        % centralities
        degree{n} = metrics.degreeCentrality{n};
        betweeness{n} = metrics.betweenessCentrality{n};
        closeness{n} = metrics.closenessCentrality{n};

        % max deg
        maxdeg{n} = max( degree{n} );                    
        imaxdeg{n} = find( degree{n} == maxdeg{n} );     
        cvcmaxdeg{n} = compVoxelCoords{n}(imaxdeg{n},:); 

        % min deg
        mindeg{n} = min( degree{n} );    
        imindeg{n} = find( degree{n} == mindeg{n} );    
        cvcmindeg{n} = compVoxelCoords{n}(imindeg{n},:);

        % max bet
        maxbet{n} = max( betweeness{n} );
        imaxbet{n} = find( betweeness{n} == maxbet{n} );
        cvcmaxbet{n} = compVoxelCoords{n}(imaxbet{n},:);

        % min bet
        minbet{n} = min( betweeness{n} );    
        iminbet{n} = find( betweeness{n} == minbet{n} );    
        cvcminbet{n} = compVoxelCoords{n}(iminbet{n},:);

        % max clo
        maxclo{n} = max( closeness{n} );
        imaxclo{n} = find( closeness{n} == maxclo{n} );    
        cvcmaxclo{n} = compVoxelCoords{n}(imaxclo{n},:);
        degmaxclo{n} = degree{n}(imaxclo{n});

        % min clo
        minclo{n} = min( closeness{n} );
        iminclo{n} = find( closeness{n} == minclo{n} );
        cvcminclo{n} = compVoxelCoords{n}(iminclo{n},:);
        degminclo{n} = degree{n}(iminclo{n});

        % farthest voxels
        [D2CFar,CVCFar,ilims,jlims,klims] = getDists2Point(compVoxelCoords{n},cvcmaxclo,1);    
        compD2CFar{n} = D2CFar;
        compCVCFar{n} = CVCFar;
        compILims{n} = ilims;
        compJLims{n} = jlims;
        compKLims{n} = klims;

        % plots    
        %plotMetricField(compVoxelCoords{n},closeness{n},'clo',n,drtVal(d),CVCFar);        
        %plotMetricField(compVoxelCoords{n},betweeness{n},'bet',n,drtVal(d));
        %plotMetricField(compVoxelCoords{n},degree{n},'deg',n,drtVal(d));              

        %-- centralities
        hdr = {'x,';'y,';'z,';'betw,';'clo,';'deg,';'kn'}; 
        hdr = hdr';        
        A = [ compVoxelCoords{n}(:,1),...
              compVoxelCoords{n}(:,2),...
              compVoxelCoords{n}(:,3),...
              betweeness{n},...
              closeness{n},...
              degree{n},...
              KN(compVoxelInds{n}) ]; 

        % writing to file
        fname = strcat('DRT_',num2str(drtVal(dv)),'Field_Cluster_',num2str(n),'_Metrics');                
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(cld,fname,'.csv'),A,'-append','delimiter',',','precision','%g');                            

        %-- max/min closeness point    

        % transform to CMG coordinates (get first in the list)
        cvcmaxcloloc = global2LocalCluster(cvcmaxclo{n}(1,1),...
                                           cvcmaxclo{n}(1,2),...
                                           cvcmaxclo{n}(1,3),...
                                           ilims(1),jlims(1),klims(1));                                               

        cvcmincloloc = global2LocalCluster(cvcminclo{n}(1,1),...
                                           cvcminclo{n}(1,2),...
                                           cvcminclo{n}(1,3),...
                                           ilims(1),jlims(1),klims(1));    

        % writing to file    
        hdr = {'Mclo,';'Mx,';'My,';'Mz,';'MxL,';'MyL,';'MzL,';...
               'mclo,';'mx,';'my,';'mz,';'mxL,';'myL,';'mzL,';...
               'im,';'iM,';'jm,';'jM,';'km,';'kM,';'s,';'R2,';'performance'}; 
        hdr = hdr';                                               
        A = [maxclo{n}(1)         ,...
             cvcmaxclo{n}(1,:)    ,...
             cvcmaxcloloc         ,...
             minclo{n}(1)         ,...
             cvcminclo{n}(1,:)    ,...
             cvcmincloloc         ,...
             ilims                ,...
             jlims                ,...
             klims                ,...
             linregr.slope{n}     ,...
             linregr.Pearson{n}   ,...
             linregr.performance{n}];

        fname = strcat('DRT_',num2str(drtVal(dv)),'Field_Cluster_',num2str(n),'_Metrics_MinMax');                                               
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');
        dlmwrite(strcat(cld,fname,'.csv'),A,'-append','delimiter',',','precision','%g');

        % farthest points    
        hdr = {'dist-L1Norm,';'xD,';'yD,';'zD,';'xDL,';'yDL,';'zDL'}; 
        hdr = hdr';                                            
        CVCFarLoc = global2LocalCluster(CVCFar(:,1),...
                                        CVCFar(:,2),...
                                        CVCFar(:,3),...
                                        ilims(1),jlims(1),klims(1));                                                                        
        A = [D2CFar,...
             CVCFar,...
             CVCFarLoc];                            

        % writing to file
        fname = strcat('DRT_',num2str(drtVal(dv)),'Field_Cluster_',num2str(n),'_FarthestPoints');                                                           
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');    
        dlmwrite(strcat(cld,fname,'.csv'),A,'-append','delimiter',',','precision','%g');

        % log data
        hdr = {'LogPhiz,';'LogRQI'}; 
        hdr = hdr';    
        A = [ linregr.logPHIZ{n},...
              linregr.logRQI{n} ];

        % writing to file
        fname = strcat('DRT_',num2str(drtVal(dv)),'Field_Cluster_',num2str(n),'_LogData');                                                           
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'');                                
        dlmwrite(strcat(cld,fname,'.csv'),A,'-append','delimiter',',','precision','%g');                
        
        % export cluster image sequence
        %cluster2image(DRT,drtVal(d),ilims,jlims,klims,fname);

    end
    
        % number of elements per cluster
        nce = cell(size(metrics.idComp));
        for ncl = 1:length(metrics.idComp)
            nce{ncl} = length(metrics.degreeCentrality{ncl});
        end
        
        % performance table per DRT
        hdr = {'cluster,';'nce,';'s,';'R2,';'performance'}; 
        hdr = hdr';    
        aux = [linregr.idComp'  ,...
               nce'              ,...
               linregr.slope'   ,...
               linregr.Pearson' ,...
               linregr.performance'];
        
        % writing to file
        fname = strcat('DRT_',num2str(drtVal(dv)),'Performance_Table');                                                           
        dlmwrite(strcat(cld,fname,'.csv'),hdr,'');                                
        dlmwrite(strcat(cld,fname,'.csv'),aux,'-append','delimiter',',','precision','%g');                
        
        clear n metrics linregr drtSt;
        
end
%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
