function [ VB, coords, inds ] = binarizeVoxels( V,val )

VB = V == val;
[I,J,K,inds]=logic2subind(VB); % Gibbon's code function
coords = [ I J K ];

end

