    %% mainVOIMetricsAnalyzer - analyzes cluster metrics in the VOI
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Jan 13th, 2016        
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
d.VOIgraphDataDependency;
d.VOIgraphMetricsDependency;

%% INPUT  

ic = 45; jc = 68; % well coordinates
P = 14;           % reservoir radius
drtVal = 14;      % DRT to get

nc = 4; % number of components to go through

%% LOAD FILES

% file name
wellfile = strcat( 'Well_I',num2str(ic),'_J',num2str(jc) );
dbase = strcat( '../mat/',wellfile,'/' );

% metrics data structure
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_MetricsData.mat'),'metrics' );
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_LinRegrData.mat'),'linregr' );

% DRT data structure
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_',wellfile,'.mat'),'VOISt' );

% load 
[~,~,~,~,KN,~,~,~,DRT] = loadMatFiles;

% dir  
cld = '../csv/clusterData/';
if exist(cld,'dir') ~= 7; mkdir(cld); end   

ncomp = numel(metrics.idComp);      % number of components
    
% nc = ncomp; % check error above comp 17 for (45,68), (26,120)!
%sumcomp = 0; % total number of voxels for the current DRT (nc components)
for n = 1:nc
    
    compVoxelInds{n} = VOISt.compVoxelInds{ metrics.idComp{n} };        
    compVoxelCoords{n} = VOISt.compVoxelCoords{ metrics.idComp{n} };    
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
                
    %sumcomp = sumcomp + VOISt.compNNodes{n};
    
    % farthest voxels
    [D2CFar,CVCFar,ilims,jlims,klims] = getDists2Point(compVoxelCoords{n},cvcmaxclo,1);    
    compD2CFar{n} = D2CFar;
    compCVCFar{n} = CVCFar;
    compILims{n} = ilims;
    compJLims{n} = jlims;
    compKLims{n} = klims;
     
    % plots
    
    %plotMetricField(compVoxelCoords{n},closeness{n},'clo',n,drtVal,CVCFar);        
    %plotMetricField(compVoxelCoords{n},betweeness{n},'bet',n,drtVal);
    %plotMetricField(compVoxelCoords{n},degree{n},'deg',n,drtVal);  
    
    % finding complementary portion of cluster for using with CMG modeling
    [outcvc,outcvi,fout,fin] = getOutVoxels(compVoxelCoords{n},n,drtVal,wellfile);
    
    % parsing file with Python
    pyi = setPyInterpreter('/usr/local/bin/python');    
    pycmd = sprintf('%s %s %s',pyi,'../py/conv_table_to_CMG.py',fout(1:end-4)); 
    [~,~] = system(pycmd);
    delete(fout); % not necessary to store
    
    % writing to file
    fname = strcat('DRT_',num2str(drtVal),'_Cluster_Metrics_',num2str(n));                
        
    % centralities
    hdr = {'x,';'y,';'z,';'betweeness,';'closeness,';'degree,';'kn'}; 
    hdr = hdr';    
    A = [ compVoxelCoords{n}(:,1),...
          compVoxelCoords{n}(:,2),...
          compVoxelCoords{n}(:,3),...
          betweeness{n},...
          closeness{n},...
          degree{n},...
          KN(compVoxelInds{n}) ]; 
        
    dlmwrite(strcat(cld,fname,'.csv'),hdr,'delimiter','');
    dlmwrite(strcat(cld,fname,'.csv'),A,'-append','delimiter',',','precision','%g');                            
    
    % max closeness point
    
    % transform to CMG coordinates (get first in the list)
    cvcmaxcloloc = global2LocalCluster(cvcmaxclo{n}(1,1),...
                                       cvcmaxclo{n}(1,2),...
                                       cvcmaxclo{n}(1,3),...
                                       ilims(1),jlims(1),klims(1));    
    dlmwrite(strcat(cld,fname,'.csv'),'P','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),cvcmaxclo{n}(1,:),'-append','delimiter',','); % global coords                     
    dlmwrite(strcat(cld,fname,'.csv'),cvcmaxcloloc,'-append','delimiter',',');      % local coords
    dlmwrite(strcat(cld,fname,'.csv'),maxclo{n}(1),'-append','precision','%g');                    
    
    % min closeness point
    
    % transform to CMG coordinates (get first in the list)        
    cvcmincloloc = global2LocalCluster(cvcminclo{n}(1,1),...
                                       cvcminclo{n}(1,2),...
                                       cvcminclo{n}(1,3),...
                                       ilims(1),jlims(1),klims(1));    
    
    dlmwrite(strcat(cld,fname,'.csv'),'Q','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),cvcminclo{n}(1,:),'-append','delimiter',','); % global coords                     
    dlmwrite(strcat(cld,fname,'.csv'),cvcmincloloc,'-append','delimiter',',');      % local coords
    dlmwrite(strcat(cld,fname,'.csv'),minclo{n}(1),'-append','precision','%g');     
    
    % farthest points
    
    CVCFarLoc = global2LocalCluster(CVCFar(:,1),...
                                    CVCFar(:,2),...
                                    CVCFar(:,3),...
                                    ilims(1),jlims(1),klims(1));    
    
    dlmwrite(strcat(cld,fname,'.csv'),'X','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),CVCFar,'-append','delimiter',',','precision','%d');
    dlmwrite(strcat(cld,fname,'.csv'),CVCFarLoc,'-append','roffset',1,'delimiter',',','precision','%d');  
          
    % distances
    dlmwrite(strcat(cld,fname,'.csv'),'D','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),D2CFar,'-append','delimiter',',','precision','%d');                    
    
    % limits
    dlmwrite(strcat(cld,fname,'.csv'),'I','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),ilims,'-append','delimiter',',','precision','%d');                    
    dlmwrite(strcat(cld,fname,'.csv'),'J','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),jlims,'-append','delimiter',',','precision','%d');                    
    dlmwrite(strcat(cld,fname,'.csv'),'K','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat(cld,fname,'.csv'),klims,'-append','delimiter',',','precision','%d');                    
    
    % Pearson, slope
    hdr = {'s','R'};
    a = [ linregr.slope{n},linregr.Pearson{n} ];
    dlmwrite(strcat(cld,fname,'.csv'),hdr,'-append','roffset',1,'delimiter',',');
    dlmwrite(strcat(cld,fname,'.csv'),a,'-append','delimiter',',','precision','%g');                        
            
    % header
    hdr = {'LogPhiz,';'LogRQI'}; 
    hdr = hdr';    
    a = [ linregr.logPHIZ{n}, linregr.logRQI{n} ];
    dlmwrite(strcat(cld,fname,'LogData.csv'),hdr,'');                                
    dlmwrite(strcat(cld,fname,'LogData.csv'),a,'-append','delimiter',',','precision','%g');                
    
    % export cluster image sequence
    cluster2image(DRT,drtVal,ilims,jlims,klims,fname);
    
    
end

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;
