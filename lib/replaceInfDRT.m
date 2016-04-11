function [ DRT ] = replaceInfDRT( DRTmatfile )
% REMOVEINFDRT replace -inf value out from DRT 3D array by 0.
%              
%
%   input: 
%       DRTmatfile: path to .mat file of DRT 3D array
%
%   output:
%       DRT: modified array
%
%  Remark: DRT = 0 indicates a cell with null porosity which 
%          may come from the FZI computation.

aux = load(DRTmatfile);
DRT = aux.DRT;             % assumes that 'DRT' is a field
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;
disp('----> -inf values switched to 0 in DRT 3D array...');

end

