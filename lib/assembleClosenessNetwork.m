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


% plotting consecutive network
scatter3(maxcloVoxelCoordsOrd(:,1),maxcloVoxelCoordsOrd(:,2),maxcloVoxelCoordsOrd(:,3),100,1:length(clID),'fill'); 
for i = 1:length(clID)-1    
        if (maxcloVoxelCoordsOrd(i,3) == maxcloVoxelCoordsOrd(i+1,3))
            cf = 'r';
        else
            cf = 'b';
        end
            hold on
            line([maxcloVoxelCoordsOrd(i,1),maxcloVoxelCoordsOrd(i+1,1)],...
                 [maxcloVoxelCoordsOrd(i,2),maxcloVoxelCoordsOrd(i+1,2)],...
                 [maxcloVoxelCoordsOrd(i,3),maxcloVoxelCoordsOrd(i+1,3)],...
                 'LineWidth',1.0,'color',cf);
        %end
    %end
end

% field's outline

% views
vws = {'xy','xz','yz','3d'};


for i = 1:length(vws)
    % vertices and faces 
    V = [1 1 1;1 220 1; 60 220 1; 60 1 1; 1 1 85; 1 220 85; 60 220 85; 60 1 85];
    F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 4 3 7 8; 2 3 7 6; 1 4 8 5];
    patch('Faces',F,'Vertices',V,'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.0);

    set(gca,'XDir','reverse','ZDir','reverse','color','none');     
    title(strcat('Max closeness network; DRT = ',num2str(drtVal)))
    xlabel('I'), ylabel('J'), zlabel('K')    
    grid off, axis equal
    colorbar
    switch vws{i}
        case 'xy'
            view([90,90]);
        case 'xz'
            view([0,0]);
        case 'yz'
            view([90,0]);
        case '3d'
            view([40,20]);
    end
    
    print('-dpdf','-r0',strcat('../figs/','Network_DRT_',num2str(drtVal),'_View_',vws{i}));

end
