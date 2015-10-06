function plotVoxelNeigh3D(ic,jc,kc,DRTV,alpha)
%{
    input: 
        ic,jc,kc: central voxel coordinate
            DRTV: 3D array with DRT values for each voxel making up the
            neighborhood
           alpha: face color transparency

    See documentation in function 'getVoxelNeigh' to know the 
    data arrangement.

    required: 
        Function 'indPatch' from Gibbon code.

%}

fig_color='w'; fig_colordef='white';
cMap=jet(250);

figure 
indPatch = 1:numel(DRTV); % getting linear indices
[F,V,C]=ind2patch(indPatch,DRTV,'v');
figuremax(fig_color,fig_colordef);
title( strcat('Voxel neighbourhood V(',num2str(ic),',',num2str(jc),',',num2str(kc),') - DRT values') );
xlabel( strcat('J=[',num2str(jc-1),',',num2str(jc),',',num2str(jc+1),']') );
ylabel( strcat('I=[',num2str(ic-1),',',num2str(ic),',',num2str(ic+1),']') );
zlabel( strcat('K=[',num2str(kc-1),',',num2str(kc),',',num2str(kc+1),']') );
hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',alpha);
axis equal; view(3); axis tight; axis vis3d; grid off;
set(gca,'XTick',[],'YTick',[],'ZTick',[])
colormap(cMap); caxis( [ min(DRTV(:)), max(DRTV(:)) ] );
colorbar

% print to file    
print('-depsc2','-r0',fullfile( '../figs/', ...
    strcat('VoxelNeigh_IC_',num2str( ic ),'_JC', ...
    num2str( jc ),'_KC_',num2str( kc ) ) ) );  

end