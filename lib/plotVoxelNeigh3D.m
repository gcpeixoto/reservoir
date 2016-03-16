function plotVoxelNeigh3D(ic,jc,kc,P,PV,alpha,ZL,flag,tname,p,varargin)
%  PLOTVOXELNEIGH3D plots the 26-neighborhood
%     input: 
%         ic,jc,kc: central voxel coordinate
%                P: ring radius  
%               PV: 3D array with property values 
%                 for each voxel making up the neighborhood
%               ZL: indices of PV whose entries are zero
%            alpha: face color transparency   
%             flag: plot only voxels with zero entry?
%                   'nonzero': yes; (default: no) 
%            tname: character title name
%                 p: print flag (true, false)
%                     
% 
%     required: 
%         Function 'ind2Patch' from Gibbon code.
% 
%  Gustavo Peixoto



if nargin < 6
    error('Missing arguments.');
elseif nargin >= 6 && nargin < 9
    ZL = [];
    flag = '';
    tname = '';
    p = false;
end
    
fig_color='w'; fig_colordef='white';
cMap=jet(250);

if strcmp(flag,'nonzero') % plot only nonzero    
    indPatch = 1:numel(PV); % getting linear indices
    indPatch = setdiff(indPatch,ZL); % only nonzero    
else % plot all
    indPatch = 1:numel(PV); 
end
    
[F,V,C]=ind2patch(indPatch,PV,'v');
figuremax(fig_color,fig_colordef);
title( strcat('Voxel neighbourhood V(',num2str(ic),',',num2str(jc),',',num2str(kc),')',tname,' - values') );

iv = num2str(ic-P:ic+P);
jv = num2str(jc-P:jc+P);
kv = num2str(kc-P:kc+P);
xlabel( strcat('J=[',jv,']') );
ylabel( strcat('I=[',iv,']') );
zlabel( strcat('K=[',kv,']') );

hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',alpha);
axis equal; view(3); axis tight; axis vis3d; grid off;
set(gca,'XTick',[],'YTick',[],'ZTick',[])
colormap(cMap); caxis( [ min(PV(:)), max(PV(:)) ] );
colorbar

if p
    % print to file    
    print('-depsc2','-r0',fullfile( '../figs/', ...
        strcat('VoxelNeigh_IC_',num2str( ic ),'_JC', ...
        num2str( jc ),'_KC_',num2str( kc ) ) ) );  
end

end