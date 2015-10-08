function [nlayers,idlayer] = getVoxelNeighNLayersByValue( VN,val )

[A,B,C] = size(VN);

MM = VN == val;

ID = [];
for a = 1:A
    for b = 1:B
        for c = 1:C
            if MM(a,b,c) == 1
                ID = [ ID; [ a b c ] ];            
            end
        end
    end
end

nlayers = length( sort( unique( ID(:,3) ) ) );
idlayer = 1:nlayers;

end

