function [icc,ncc,mncc,isol,nisol] = plotCompHist( val, cvc, lim, w, pad, varargin )
%   PLOTCOMPHIST plots the histogram of distribution of the connected
%                components of a DRT network
%   input:
%           val: DRT value;
%           cvc: cell containing the component's voxel coordinates
%           lim: number of components to plot in the histogram
%             w: width of the histogram bar
%           pad: padding vector for the plot: [ padx, pady ]
%   output:
%           icc: index of the connected component
%           ncc: number of voxels in that component
%          mncc: mean number of voxels for that network
%          isol: index of isolated voxels (components with 1 voxel only)
%         nisol: number of isolated voxels
%             h: graphic handle
%
%

if ~isa(cvc,'cell');
    error('Argument must be a cell.');
end

nComps = size(cvc,2); % number of components

if nargin == 2    
    lim = nComps;
    w = 0.8;
    pad = [0.5,5];
elseif nargin == 3
    w = 0.8;
    pad = [0.5,5];
elseif nargin == 4
    pad = [0.5,5];
end

% table for bar plot
icc = zeros(nComps,1);
ncc = 0*icc;
for i = 1:nComps
    icc(i) = i;
    ncc(i) = length( cvc{i} );
end

% isolated components
nisol = 0; 
isol = [];
idisol = find(ncc == 1);
if ~isempty(idisol)
    isol = icc(idisol);
    nisol = length(isol);
end

icc = icc(1:lim); ncc = ncc(1:lim); % reducing
mncc = round(mean(ncc));   

% bar plot
figure
bar(icc,ncc,w,'k'); 
%%---- plot line of mean value per component
%hold on
%line([ 0-padx,max(icc)+padx ], [ mncc,mncc ],'Color', [1 0 0],'LineWidth',2); % maximum line
% TODO improve labels if necessary
str = strcat('DRT=',num2str(val) );
legend({str});
xlabel('Connected components');
ylabel('Voxel numbers');
set(gca,'TickDir','out');
axis([0-pad(1) max(icc)+pad(1) 0-pad(2) max(ncc)+pad(2)])

end

