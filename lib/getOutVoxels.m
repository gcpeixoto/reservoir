function [ outcvc, outcvi, fout, fin ] = getOutVoxels( cvc,comp,drtVal,wellname,ext,varargin )
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
%        fout: output file
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
m = min(cvc,[],1); 
M = max(cvc,[],1);
[im,jm,km,iM,jM,kM] = deal(m(1),m(2),m(3),M(1),M(2),M(3));

% [x y z] of all voxels
[I,J,K] = ndgrid(im:iM, jm:jM, km:kM);
allcvc = [I(:),J(:),K(:)];

% operation of exclusion: out = all - in
[outcvc,outcvi] = setdiff(allcvc,cvc,'rows');
           
% write to file
dn = strcat('../txt/',wellname);
if exist(dn,'dir') ~= 7
    mkdir(dn);
end

if exist('../csv/clusterData/','dir') ~= 7
    mkdir('../csv/clusterData/');
end

fout = strcat( dn,'/','OutCoords_DRT_',num2str(drtVal),... 
                   '_Cluster_',num2str(comp),'_',wellname,ext );    

fin = strcat( '../csv/clusterData/','InCoords_DRT_',num2str(drtVal),... 
                   '_Cluster_',num2str(comp),'_',wellname,'.csv' );    
               
               
dlmwrite(fout,['x' 'y' 'z'],'delimiter',' ');
dlmwrite(fout,outcvc,'-append','delimiter',' '); 

dlmwrite(fin,['x' 'y' 'z'],'delimiter',',');
dlmwrite(fin,cvc,'-append','delimiter',','); 

end

