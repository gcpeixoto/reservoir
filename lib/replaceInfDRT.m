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

if ischar(DRTmatfile)   % path to .mat file
    aux = load(DRTmatfile);
    DRT = aux.DRT;             % assumes that 'DRT' is a field
else                    % numeric
    DRT = DRTmatfile;
end

DRT(DRT == -Inf) = 0.0; % eliminating -Inf

    disp('----> -inf values switched to 0 in DRT 3D array...');

end

