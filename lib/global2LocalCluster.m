function loc = global2LocalCluster( xglob,yglob,zglob,xm,ym,zm)
% GLOBAL2LOCALCLUSTER coordinate transform for cluster coordinates
%
%   input:
%       xglob,yglob,zglob: voxel's global coordinates (whole field)
%                xm,ym,zm: cluster's initial limits                        
%   output:
%       loc: voxel local coordinates at cluster (1x3)
%

xloc = xglob - xm + 1;
yloc = yglob - ym + 1;
zloc = zglob - zm + 1;

loc = [xloc,yloc,zloc];

end

