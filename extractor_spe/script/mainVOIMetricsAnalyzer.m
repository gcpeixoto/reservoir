%% mainVOIMetricsAnalyzer
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Jan 13th, 2016        
%             
%   description: finds points of minimum and maximum centralities, plots 
%                the fields and export data to files

clear all; close all; clc; format long;
set(0,'DefaultAxesFontSize',18);  
%% 

% well
ic = 26; jc = 120;

wellfile = strcat( 'Well_I',num2str(ic),'_J',num2str(jc) );
dbase = strcat( '../mat/',wellfile,'/' );

drtVal = '16'; % DRT to get

% metrics data structure
aux = load( strcat(dbase,'VOI_DRT_',drtVal,'_MetricsData.mat') );
metrics = aux.metrics;

aux = load( strcat(dbase,'VOI_DRT_',drtVal,'_LinRegrData.mat') );
linregr = aux.linregr;

% DRT data structure
aux = load( strcat(dbase,'VOI_DRT_',drtVal,'_',wellfile,'.mat') );
VOISt = aux.VOISt;
val = VOISt.value;

ncomp = numel(metrics.idComp);      % number of components

nc = 5;
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
       
    plotMetricField(compVoxelCoords{n},closeness{n},'clo',n,val,CVCFar);        
    %plotMetricField(compVoxelCoords{n},betweeness{n},'bet',n,val);
    %plotMetricField(compVoxelCoords{n},degree{n},'deg',n,val);  
    
    
    % writing to file
    fname = strcat('DRT_',num2str(val),'_Info_Cluster_',num2str(n));            
            
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
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'P','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcmaxclo{n}(1,:),'-append','delimiter',',');                    
    
    % min closeness point
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'Q','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),cvcminclo{n}(1,:),'-append','delimiter',',');
    
    % farthest points
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),'X','-append','roffset',1,'delimiter',',');                    
    dlmwrite(strcat('../csv/subdomain/',fname,'.csv'),CVCFar,'-append','delimiter',',','precision','%d');                    
          
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
