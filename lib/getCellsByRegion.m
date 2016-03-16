function [ cvcg,cvig,cvco,cvio,cvcw,cviw ] = getCellsByRegion( cvc,zgo,zow )
%GETCELLSBYREGION locates cells inside regions of gas, oil and water
%
%   input:
%        cvc: cluster's global cell coords
%        zgo: z coord of gas-oil contact
%        zow: z coord of oil-water contact
%   output:
%       cvcg: cluster's global cell coords (gas region)
%       cvig: cluster's local cell coords (gas region)

% locating cells inside each region
ig =  cvc(:,3) <= zgo ;
io =  cvc(:,3) > zgo & cvc(:,3) <= zow ;
iw =  cvc(:,3) > zow ;

% global coords
cvcg = cvc(ig,:);
cvco = cvc(io,:);
cvcw = cvc(iw,:);

% cluster's local indices
cvig = find(ig); 
cvio = find(io);
cviw = find(iw);

end

