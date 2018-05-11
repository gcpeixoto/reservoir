function dv = setNeighDist(neigh)
%SETNEIGH chooses the neighbourhood connectivity
%
% input: 
%       neigh = connectivity structuring element ('6','26')
%
% output: 
%          dv = distance to use as criterion for connectivity
%

% ischar 
assert(ischar(neigh),'neigh must be a char: %s or %s','6','26');

switch neigh
    
    % 6 neighbours (3D von Neumman's neighborhood) => distance <=1
    %       _
    %     _|_|_    |
    %    |_|_|_| - o -
    %      |_|     |
    %       
    %   i-1 i i+1
    %
    case '6'       
        dv = 1;        
    
        
    % 26 neighbours (3D Moore's neighborhood) => distance <= sqrt(3)
    %     _ _ _
    %    |_|_|_|  \ | /
    %    |_|_|_|  - o -
    %    |_|_|_|  / | \
    %       
    %   i-1 i i+1
    %
    case '26'
        dv = sqrt(3);
        
end
        

end

