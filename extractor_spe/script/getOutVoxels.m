function [ outcvc, outcvi ] = getOutVoxels( cvc,comp,drtVal,wellname,ext,varargin )
% GETOUTVOXELS extract the complementary set of a cluster bounded by 
%              a box
%
%   input:
%         cvc: cluster's voxel coordinates (array nx3)
%        comp: component number (int)
%      drtVal: DRT value
%         ext: file extension (default: .txt)
%   output:
%      outcvc: complementary of cluster cluster's coordinates
%      outcvi: voxel indices in relation to box's indices
%
%   Remark: .txt file required to be imported into CMG in order to 
%           characterize the cluster's irregular geometry
%
%
% CLUSTER'S SURROUNDINGS
% ======================
%              _____________________
%             | |_|           o     |???? cluster bounding box
%             | |_|  o       o      |
%             | |_|   o    o     o  |
%             | |_|_ _   o o        |
%             | |_|_|i|_          o | 
%             |     |i|_|_ _ _ _ _ _|
%             |  o  |i| |_|_|_|_|_|_|   i: in  (cluster voxels)
%             | o   |i|  o     o    |   o: out (complementary voxels)
%             |_____|_|_____________|   out = all - in
%
% 

if nargin == 4        
    ext = '.txt';         
end

% cluster bounding box limits
im=min(cvc(:,1)); iM=max(cvc(:,1)); 
jm=min(cvc(:,2)); jM=max(cvc(:,2)); 
km=min(cvc(:,3)); kM=max(cvc(:,3));

% cuboid grid formed by all voxels in the cluster's bounding box
I = im:iM;
J = jm:jM;
K = km:kM;

% [x y z] of all voxels
allcvc = [];
for i = 1:length(I)
    for j = 1:length(J)
        for k = 1:length(K)            
            allcvc = [ allcvc; [I(i),J(j),K(k)] ];  
        end
    end
end

% operation of exclusion: out = all - in
[outcvc,outcvi] = setdiff(allcvc,cvc,'rows');
           
% write to file
fout = strcat( '../txt/','OutCoords_DRT_',num2str(drtVal),... 
                   '_Cluster_',num2str(comp),'_',wellname,ext );    

dlmwrite(fout,['x' 'y' 'z'],'delimiter',' ');
dlmwrite(fout,outcvc,'-append','delimiter',' '); 

end

