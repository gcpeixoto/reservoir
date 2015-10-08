function plotVoxelNeigh2DLayer( ic,jc,layer,DRTV)
%{
    input: 
        ic,jc,kc: central voxel coordinate
            DRTV: 3D array with DRT values for each voxel making up the
            neighborhood
           layer: neighborhood's layer             
    

    Gustavo Peixoto, Oct 8 2015, @UFPB
%}

if layer < 1 || layer > size(DRTV,3)
    error('Layer value out of range.'); 
end

cMap=jet(250);

figure
L_plot = false(size(DRTV)); 
L_plot(:,:,layer) = 1;
ind = find(L_plot);
[F,V,C]=ind2patch(ind,DRTV,'sk'); % creating patch data for layer
patch('Faces',F,'Vertices',V,'EdgeColor','k', 'CData',C,'FaceColor','flat','FaceAlpha',1);

title( strcat('Voxel neighbourhood slice: V(:,:,','k=',num2str(layer),')' ) );
xlabel( strcat('J=[',num2str(jc-1),',',num2str(jc),',',num2str(jc+1),']') );
ylabel( strcat('I=[',num2str(ic-1),',',num2str(ic),',',num2str(ic+1),']') );
set(gca,'XTick',[],'YTick',[])
colormap(cMap); caxis( [ min(DRTV(:)), max(DRTV(:)) ] );
colorbar 

% print to file    
print('-depsc2','-r0',fullfile( '../figs/', ...
    strcat('VoxelNeigh_IC_',num2str( ic ),'_JC', ...
    num2str( jc ),'_KSlice_',num2str( layer ) ) ) );  

end

