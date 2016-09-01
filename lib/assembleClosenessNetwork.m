function [maxcloVoxelCoordsOrd,origID,clID,maxclo] = assembleClosenessNetwork(drtVal,dirpath)
%     ASSEMBLECLOSENESSNETWORK assembles the 3D network of high-performance
%                              (HP) clusters per DRT based on closeness
%                              centrality.
% 
%     input:
%         drtVal: DRT value         
%         dirpath: directory path to .mat files. (default: '../mat/Field/') 
% 
%     output:
%         maxcloVoxelCoordsOrd: array with 3D grid coordinates of the
%                               max closeness cell of each cluster (nc x 3)
%         origID: array with original HP clusters IDs in the order modified 
%                 after calling 'sortrow' to sort the coordinates in terms
%                 of z-depth (lower to higher) (nc x 1)
%         clID: HP cluster IDs in the originial order of search (nc x 1)
%         maxclo: array with max closeness per HP cluster (nc x 1)
% 
%     - REMARK
%     
%     This strategy is, at first glance, ONLY applicable to 
%     .dat (.txt) files that follow the output format of CMG software. 
%     The current version is tested for the SPE Project 2 
%     data arrangement. It will might be useful for other 
%     files as well, but further investigation is required. 
% 
%     What is done here: the original files are 'unrolled' 
%     and rearranged into 3D arrays (I,J,K), thus allowing a 
%     better accessibility of the information. In the end, the reservoir
%     is a 3D structure with P(I,J,K) evaluated 
%     in the fourth dimension. 
% 
%     RESERVOIR SCHEME
%     ================
% 
%           K
%          /
%         /____________
%        /|            |
%       /_|__________  |
%      /|            | |
%     /_|__________  |_|      <------ (i,j,k = K): bottom layer
%     |            | |
%     |            |_|        <------ (i,j,k = 2): depth layer
%     |            | 
%     |____________|___ J     <------ (i,j,k = 1): reservoir's surface
%     |    
%     |
%      I
%
% -----------------------
% Dr. Gustavo Peixoto


if isempty(dirpath), dirpath = '../mat/Field/'; end

% load structures into workspace and store
load( strcat(dirpath,'DRT_',num2str(drtVal),'.mat'),'drtSt' );
load( strcat(dirpath,'DRT_',num2str(drtVal),'_MetricsData.mat'),'metrics' );
load( strcat(dirpath,'DRT_',num2str(drtVal),'_LinRegrData.mat'),'linregr' );        

% seeks only high-performance clusters 
clID = [];     
for n = 1:length(metrics.idComp)        
    if linregr.performance{n} == 1 
        clID(n) = metrics.idComp{n};                                    
    else
        clID(n) = 0;
    end            
end

% purges low-performance and checks if remains someone
clID = clID(clID~=0); 
assert(~isempty(clID),sprintf(['No HP cluster found for DRT = %d. '...
                                'Network cannot be formed.'],drtVal));

% getting max closeness point
maxclo = cell(1,length(clID));
for cid = 1:length(clID)
    maxclo{cid} = metrics.closenessCentrality{clID(cid)}; % separating output
    maxcloVoxelCoords{cid} = metrics.maxClosenessVoxelCoords{clID(cid)}(1,:); % coords
end
maxclo=cellfun(@max,maxclo); % max closeness 

%{
    gets the max closeness point coords and sort them according to the 
    z value (3rd column).
    nodeID is the original cluster numbering obtained by sweeping the
    'metrics' struct.
%}
maxcloVoxelCoords = reshape(cell2mat(maxcloVoxelCoords)',[3,numel(maxcloVoxelCoords)])';
[maxcloVoxelCoordsOrd,origID] = sortrows(maxcloVoxelCoords,3);

%{

scatter3(maxcloVoxelCoordsOrd(:,2),maxcloVoxelCoordsOrd(:,1),maxcloVoxelCoordsOrd(:,3),100,maxcloVoxelCoordsOrd(:,3),'fill'); 
for i = 1:size(maxcloVoxelCoordsOrd,1)
    for j = i:size(maxcloVoxelCoordsOrd,1)
        if (i~=j)
            hold on
            line([maxcloVoxelCoordsOrd(i,2),maxcloVoxelCoordsOrd(j,2)],...
                 [maxcloVoxelCoordsOrd(i,1),maxcloVoxelCoordsOrd(j,1)],...
                 [maxcloVoxelCoordsOrd(i,3),maxcloVoxelCoordsOrd(j,3)],...
                 'LineWidth',0.1,'color','k');
        end
    end
end
set(gca,'ZDir','reverse'); view([10,13])
print('-dpdf','-r0',fullfile( '../figs/network-all'))

%}

end
