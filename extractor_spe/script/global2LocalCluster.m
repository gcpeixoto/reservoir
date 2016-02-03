function loc = global2LocalCluster( xglob,yglob,zglob,zglobsup,ic,jc,P )
% GLOBAL2LOCALCLUSTER coordianate transform for cluster coordinates
%
%   input:
%       xglob,yglob,zglob: voxel's global coordinates (whole field)
%                zglobsup: global z coordinate of cluster's uppermost voxel
%                       P: tier's radius (reservoir)
%   output:
%       loc: voxel local coordinates at cluster (1x3)
%
% This method outputs the depth local coordinates in terms of the cluster
% height. In CMG, zglobsup (field) renders the z = 1 (local). Then, we 
% need work with a grid of coordinates based on the CLUSTER'S DOMAIN.
%
% REMARK: method valid for the symmetric reservoirs that are opened
%         from the central well's coordinates (ic,jc).
%        
%
%
%              -P         (ic,jc)         +P
%               <------------|------------>
%               ___________________________
%              |                           |
%              |                           |
%              |   _ _ _   __ zglobsup     |  _ zglobsup => zloc = 1
%              |  |_|_|_|                  |
%              |    |_|_|  __ zglob        |
%              |      |_|  __ zglobinf     |  _ zglobinf => zloc = zcluster              
%              |                           |
%              |___________________________| 
%
% transform law

xloc = P+1 - (ic - xglob);
yloc = P+1 - (jc - yglob);
zloc = zglob - zglobsup + 1;

loc = [xloc,yloc,zloc];

end

