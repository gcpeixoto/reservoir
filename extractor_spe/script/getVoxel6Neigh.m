function [ vc6n, vi6n, phi6n, kn6n ] = getVoxel6Neigh( vc, vList, viList, PHI, KN )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

if numel(vc) ~= 3 || size(vList,2) ~= 3 || size(viList,2) ~= 1
    error('Check arguments.')
end

vc6n = zeros(6,3);
vi6n = zeros(6,1);
phi6n = zeros(6,1);
kn6n = zeros(6,1);

for  i = 1:length(vList)
    
    % +x
    aux = [ vc(1)+1, vc(2), vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(1,:) = aux;        
        vi6n(1) = ind;
        phi6n(1) = PHI(ind);
        kn6n(1) = KN(ind);
    end
    
    % -x
    aux = [ vc(1)-1, vc(2), vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(2,:) = aux;        
        vi6n(2) = ind;
        phi6n(2) = PHI(ind);
        kn6n(2) = KN(ind);
    end
    
    % +y
    aux = [ vc(1), vc(2)+1, vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(3,:) = aux;        
        vi6n(3) = ind;
        phi6n(3) = PHI(ind);
        kn6n(3) = KN(ind);
    end
    
    % -y
    aux = [ vc(1), vc(2)-1, vc(3) ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(4,:) = aux;        
        vi6n(4) = ind;
        phi6n(4) = PHI(ind);
        kn6n(4) = KN(ind);
    end
    
    % +z
    aux = [ vc(1), vc(2), vc(3)+1 ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(5,:) = aux;        
        vi6n(5) = ind;
        phi6n(5) = PHI(ind);
        kn6n(5) = KN(ind);
    end
    
    % -z
    aux = [ vc(1), vc(2), vc(3)-1 ];
    ind = strmatch( aux, vList );
    if ~isempty(ind)
        vc6n(6,:) = aux;        
        vi6n(6) = ind;
        phi6n(6) = PHI(ind);
        kn6n(6) = KN(ind);
    end
    
end
    

