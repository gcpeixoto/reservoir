function plotVoxelGraphComp( DRTMat, inds, val, idComp )
%{
    PLOTVOXELGRAPHCOMP - plot the voxels with
                   DRT=val for all the reservoir
    input: 
        DRTMat: 3D array of DRT values
          inds: indices of the voxels that form 
                a graph component (isolated or connected)
           val: DRT value 
        idComp: indice of the graph component (integer that refer to 
                a unique isolated node or a family of connected nodes)
%}


figure 
[F,V,C]=ind2patch(inds,DRTMat,'v');
        
hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.8);
axis equal; view(3); axis tight; axis vis3d; grid off;        
title( strcat( 'Voxel Component: ',num2str(idComp),'- DRT ', num2str( val ) ) );
xlabel('J');
ylabel('I');
zlabel('K');

%print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
%        strcat('Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ) ) ) );  
    
%saveas(gcf, strcat( '../figs/graphPath/', ...
%        'Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ),'.fig' ) );  

end

