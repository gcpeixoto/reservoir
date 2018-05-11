%% mainFieldExportCluster2CMG - export table of deactivated out cluster cells
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Aug 31th, 2016        
%             
%   description: write output files having the Cartesian coordinates 
%                of deactivated cells to be used with CMG IMEX.
%                Since the cluster's domain is a irregular geometry,
%                in order to have it, it is necessary to deactivate
%                all the cells that do not make up this geometry.
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
d.fieldGraphMetricsDependency;

%% INPUT

%{
    How to use:

        - drtVal: array of DRT values to get. Unique or list.                
        - comps: array of cluster IDs.                                      

    REMARK: if the user choose, e.g. drtVal = [1,4] and comps = [3,10,15],
    the code will search for the clusters 3,10,15 for DRT = 1,4. However, 
    the user must verify if such clusters exist for the respective DRT. 
    
    Furthermore, the utility of this method relies on selecting only special
    clusters for study, like the high-performance ones. Therefore, it is
    recommended that the user selects only one DRT value and the desired 
    cluster IDs. 

    This method was tailored for cluster analysis over the whole FIELD 
    of SPE10 dataset.
                        
%}

drtVal = 14;
comps = [2,3,8];

load('../mat/DRT.mat','DRT');

for dv = 1:length(drtVal)

    % metrics data structure
    load( strcat('../mat/Field/DRT_',num2str(drtVal(dv)),'_MetricsData.mat'),'metrics' );
    load( strcat('../mat/Field/DRT_',num2str(drtVal(dv)),'_LinRegrData.mat'),'linregr' );

    % DRT data structure
    load( strcat('../mat/Field/DRT_',num2str(drtVal(dv)),'.mat'),'drtSt' );

    for dn = 1:length(comps)% cluster loop    
        
        %{
        % plotting
        cvi{dn} = drtSt.compVoxelInds{ metrics.idComp{comps(dn)} };         
        plotVoxelGraphComp(DRT,cvi{dn},drtVal(dv),comps(dn),1.0,'default','p',[-37.5,30],'pdf');
        clf;
        %}
                
        compVoxelCoords{dn} = drtSt.compVoxelCoords{ metrics.idComp{comps(dn)} };        
        
        % finding complementary portion of cluster for using with CMG modeling
        [outcvc,outcvi,fout,fin] = getOutVoxels(compVoxelCoords{dn},...
                                                comps(dn),drtVal(dv),'Field');
                                            
        % parsing file with Python
        fprintf('--> Parsing file %s with Python...\n',fout)
        pyi = setPyInterpreter('/usr/bin/python');    
        pycmd = sprintf('%s %s %s',pyi,'../py/conv_table_to_CMG.py',fout(1:end-4)); 
        [status,~] = system(pycmd);
        assert(status == 0,'Problem found while calling the Python interpreter. Check the path.');
        
        % converting .txt -> .inc
        fprintf('--> Converting to .inc to be included into IMEX .dat file.\n')
        [path,name,ext] = fileparts(fout);
        fparsed = strcat(path,filesep,name,'_Parsed',ext);
        
        if ~exist(strcat(path,filesep,'inc'),'dir'), mkdir(path,'inc'), end;        
        fnew = strcat(path,filesep,name,'_Parsed','.inc');
        shcmd = strcat(['mv ',fparsed,' ',fnew,'; mv ',fnew,' ',path,filesep,'inc']); 
        [~,~] = system(shcmd);
                       
        delete(fout); % not necessary to store
        
    end

end