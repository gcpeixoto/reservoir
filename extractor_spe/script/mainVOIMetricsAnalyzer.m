%% mainVOIMetricsAnalyzer
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Jan 13th, 2016        
%             
%   description: finds points of minimum and maximum centralities, plots 
%                the fields and export data to files


  
%% Default

clear all; close all; clc;
activateLog(mfilename);  % log
setOptions;

%% 

ic = 45; jc = 68; % well coordinates
P = 14;           % reservoir radius
drtVal = 13;      % DRT to get

% file name
wellfile = strcat( 'Well_I',num2str(ic),'_J',num2str(jc) );
dbase = strcat( '../mat/',wellfile,'/' );

% metrics data structure
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_MetricsData.mat'),'metrics' );
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_LinRegrData.mat'),'linregr' );

% DRT data structure
load( strcat(dbase,'VOI_DRT_',num2str(drtVal),'_',wellfile,'.mat'),'VOISt' );

ncomp = numel(metrics.idComp);      % number of components

nc = 3;
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
    
    % finding complemetary portion of cluster for using with CMG modeling
    [outcvc,outcvi,fout] = getOutVoxels(compVoxelCoords{n},n,drtVal,wellfile);
    
    % parsing file with Python
    pyi = setPyInterpreter('/usr/local/bin/python');    
    pycmd = sprintf('%s %s %s',pyi,'../py/conv_table_to_CMG.py',fout(1:end-4)); 
    [~,~] = system(pycmd);
    delete(fout); % not necessary to store
    
    % writing to file
    fname = strcat('DRT_',num2str(drtVal),'_Info_Cluster_',num2str(n));                
        
    % centralities
    hdr = {'x','y','z','b','c','d'};    
    A = [ compVoxelCoords{n}(:,1),...
          compVoxelCoords{n}(:,2),...
          compVoxelCoords{n}(:,3),...
          betweeness{n},...
          closeness{n},...
          degree{n} ];                        
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),hdr,'delimiter',',');
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),A,'-append','delimiter',',','precision','%g');                        
    
    % max closeness point
    
    % transform to CMG coordinates
    cvcmaxcloloc = global2LocalCluster(cvcmaxclo{n}(1),...
                                       cvcmaxclo{n}(2),...
                                       cvcmaxclo{n}(3),...
                                     klims(1),ic,jc,P);    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'P','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcmaxclo{n}(1,:),'-append','delimiter',','); % global coords                     
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcmaxcloloc,'-append','delimiter',',');      % local coords
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),maxclo{n}(1),'-append','precision','%g');                    
    
    % min closeness point
    
    % transform to CMG coordinates
    cvcmincloloc = global2LocalCluster(cvcminclo{n}(1),...
                                       cvcminclo{n}(2),...
                                       cvcminclo{n}(3),...
                                     klims(1),ic,jc,P);    
    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'Q','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcminclo{n}(1,:),'-append','delimiter',','); % global coords                     
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcmincloloc,'-append','delimiter',',');      % local coords
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),minclo{n}(1),'-append','precision','%g');     
    
    % farthest points
    
    CVCFarLoc = global2LocalCluster(CVCFar(:,1),...
                                    CVCFar(:,2),...
                                    CVCFar(:,3),...
                              klims(1),ic,jc,P); 
    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'X','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),CVCFar,'-append','delimiter',',','precision','%d');
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),CVCFarLoc,'-append','roffset',1,'delimiter',',','precision','%d');  
          
    % distances
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'D','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),D2CFar,'-append','delimiter',',','precision','%d');                    
    
    % limits
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'I','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),ilims,'-append','delimiter',',','precision','%d');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'J','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),jlims,'-append','delimiter',',','precision','%d');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'K','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),klims,'-append','delimiter',',','precision','%d');                    
    
    % Pearson, slope
    hdr = {'s','R'};
    a = [ linregr.slope{n},linregr.Pearson{n} ];
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),hdr,'-append','roffset',1,'delimiter',',');
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),a,'-append','delimiter',',','precision','%g');                        
        
               
end
