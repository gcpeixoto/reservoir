function [ vc6n, vi6nloc, vi6nglob, f6n ] = getVoxel6Neigh( vc, vList, viList, F)
% GETVOXEL6NEIGH gets the voxel coordinates and indices in a list
%                of voxels (cluster) according to a 6-neighbourhood 
%                criterion 
%
%   input:
%          vc: seed-voxel coordinates
%       vList: cluster voxel's coordinates
%      viList: cluster voxel's indices
%           F: 3D array of the scalar field
%
%   output:
%        vc6n: 6-neighbour voxel coordinates
%     vi6nloc: 6-neighbour local voxel indices
%    vi6nglob: 6-neighbour global voxel indices
%         f6n: 6-neighbour field value
%
%
% 6-neighbour array: [+x,-x,+y,-y,+z,-z]'
%

if numel(vc) ~= 3 || size(vList,2) ~= 3 || size(viList,2) ~= 1
    error('Check arguments.')
end

vc6n = zeros(6,3);
vi6nloc = zeros(6,1);
vi6nglob = zeros(6,1);
f6n = zeros(6,1);

% sweep cluster to find neighbours  

% +x 
aux = [ vc(1)+1, vc(2), vc(3) ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(1,:) = aux;        
    vi6nloc(1) = ind;
    vi6nglob(1) = viList(ind);
    f6n(1) = F( viList(ind) );        
end

% -x
aux = [ vc(1)-1, vc(2), vc(3) ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(2,:) = aux;        
    vi6nloc(2) = ind;
    vi6nglob(2) = viList(ind);
    f6n(2) = F( viList(ind) );        
end

% +y
aux = [ vc(1), vc(2)+1, vc(3) ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(3,:) = aux;        
    vi6nloc(3) = ind;
    vi6nglob(3) = viList(ind);
    f6n(3) = F( viList(ind) );        
end

% -y
aux = [ vc(1), vc(2)-1, vc(3) ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(4,:) = aux;        
    vi6nloc(4) = ind;
    vi6nglob(4) = viList(ind);
    f6n(4) = F( viList(ind) );        
end

% +z
aux = [ vc(1), vc(2), vc(3)+1 ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(5,:) = aux;        
    vi6nloc(5) = ind;
    vi6nglob(5) = viList(ind);
    f6n(5) = F( viList(ind) );        
end

% -z
aux = [ vc(1), vc(2), vc(3)-1 ];
ind = find( ismember( vList,aux,'rows' ) );
if ~isempty(ind)
    vc6n(6,:) = aux;        
    vi6nloc(6) = ind;
    vi6nglob(6) = viList(ind);
    f6n(6) = F( viList(ind) );        
end
    
id = vi6nloc ~= 0;
vi6nloc = vi6nloc(id);
vi6nglob = vi6nglob(id);
vc6n = vc6n(id,:);
f6n = f6n(id);


end
    

