function plotDRTField( DRT, val )
%{
    PLOTDRTFIELD - plot the voxels with
                   DRT=val for all the reservoir
    input: 
        DRT: 3D array of DRT values
        val: DRT value chosen
%}


figure 
[ ~, ~, ind ] = maskVoxels( DRT,val );
[F,V,C]=ind2patch(ind,DRT,'v');
        
hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.8);
axis equal; view(3); axis tight; axis vis3d; grid off;        
title( strcat( 'Path - DRT ', num2str( val ) ) );
xlabel('J');
ylabel('I');
zlabel('K');

print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
        strcat('Voxel_Regions_DRT_',num2str( val ) ) ) );  

end

