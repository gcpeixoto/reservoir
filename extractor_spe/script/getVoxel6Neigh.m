function [ vc6n, vi6n, f6n ] = getVoxel6Neigh( vc, vList, viList, F)
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
%        vi6n: 6-neighbour voxel indices
%         f6n: 6-neighbour field value
%
%
% 6-neighbour array: [+x,-x,+y,-y,+z,-z]'
%

if numel(vc) ~= 3 || size(vList,2) ~= 3 || size(viList,2) ~= 1
    error('Check arguments.')
end

vc6n = zeros(6,3);
vi6n = zeros(6,1);
f6n = zeros(6,1);

% sweep cluster
for  i = 1:length(vList)
    
    % +x 
    aux = [ vc(1)+1, vc(2), vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(1,:) = aux;        
        vi6n(1) = ind;
        f6n(1) = F(ind);        
    end
    
    % -x
    aux = [ vc(1)-1, vc(2), vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(2,:) = aux;        
        vi6n(2) = ind;
        f6n(2) = F(ind);        
    end
    
    % +y
    aux = [ vc(1), vc(2)+1, vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(3,:) = aux;        
        vi6n(3) = ind;
        f6n(3) = F(ind);        
    end
    
    % -y
    aux = [ vc(1), vc(2)-1, vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(4,:) = aux;        
        vi6n(4) = ind;
        f6n(4) = F(ind);        
    end
    
    % +z
    aux = [ vc(1), vc(2), vc(3)+1 ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(5,:) = aux;        
        vi6n(5) = ind;
        f6n(5) = F(ind);        
    end
    
    % -z
    aux = [ vc(1), vc(2), vc(3)-1 ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(6,:) = aux;        
        vi6n(6) = ind;
        f6n(6) = F(ind);        
    end
    
end
id = vi6n ~= 0;
vi6n = vi6n(id);
vc6n = vc6n(id,:);
f6n = f6n(id);


end
    

