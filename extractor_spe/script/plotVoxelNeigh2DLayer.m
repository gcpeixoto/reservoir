function plotVoxelNeigh2DLayer( ic,jc,kc,DRTV,layer)
%{
    input: 
        ic,jc,kc: central voxel coordinate
            DRTV: 3D array with DRT values for each voxel making up the
            neighborhood
           layer: neighborhood's layer 
            ( 1 = layer kc-1; 2 = layer k; 3 = layer kc+1 )

    See documentation in function 'getVoxelNeigh' to know the 
    data arrangement.

%}

if layer == 1
    k = kc-1;
elseif layer == 2
    k = kc;
elseif layer == 3
    k = kc + 1;
else
    error('k value is 1,2 or 3.');
end

cMap=jet(250);

figure
L_plot = false(size(DRTV)); 
L_plot(:,:,layer) = 1;
ind = find(L_plot);
[F,V,C]=ind2patch(ind,DRTV,'sk'); % creating patch data for layer
patch('Faces',F,'Vertices',V,'EdgeColor','k', 'CData',C,'FaceColor','flat','FaceAlpha',1);

title( strcat('Voxel neighbourhood slice: V(:,:,','k=',num2str(k),')' ) );
xlabel( strcat('J=[',num2str(jc-1),',',num2str(jc),',',num2str(jc+1),']') );
ylabel( strcat('I=[',num2str(ic-1),',',num2str(ic),',',num2str(ic+1),']') );
set(gca,'XTick',[],'YTick',[])
colormap(cMap); caxis( [ min(DRTV(:)), max(DRTV(:)) ] );
colorbar 

% print to file    
print('-depsc2','-r0',fullfile( '../figs/', ...
    strcat('VoxelNeigh_IC_',num2str( ic ),'_JC', ...
    num2str( jc ),'_KSlice_',num2str( k ) ) ) );  

end

