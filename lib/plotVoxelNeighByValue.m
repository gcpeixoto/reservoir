function plotVoxelNeighByValue( ic,jc,kc,P,PV,val,alpha,tname,p)
% PLOTVOXELNEIGHBYVALUE
%     input: 
%         ic,jc,kc: central voxel coordinate
%                P: ring radius  
%               PV: 3D array with property values 
%                 for each voxel making up the neighborhood
%               ZL: indices of PV whose entries are zero
%            alpha: face color transparency               
%            tname: variable name to append in file name
%                p: print flag (true,false)
%                     
% 
%     required: 
%         Function 'ind2patch' from Gibbon code.
% 
%  Gustavo Peixoto




[A,B,C] = size(PV);

MM = PV == val;

% selecting indices with given value
ID = [];
ID2 = [];
for a = 1:A
    for b = 1:B
        for c = 1:C
            if MM(a,b,c) == 1
                ID = [ ID; [ a b c ] ];
            else
                ID2 = [ ID2; [ a b c ] ];
            end
        end
    end
end

% interest voxels
IND = [];
for p = 1:size(ID,1)    
    IND = [IND, sub2ind(size(PV),ID(p,1),ID(p,2),ID(p,3)) ]; % get linear indices
end

% complementary voxels
IND2 = [];
for p = 1:size(ID2,1)    
    IND2 = [IND2, sub2ind(size(PV),ID2(p,1),ID2(p,2),ID2(p,3)) ]; % get linear indices
end

fig_color='w'; fig_colordef='white';
cMap  = jet(250);

[F,V,C]=ind2patch(IND,PV,'v');
figuremax(fig_color,fig_colordef);
title( strcat('Voxel neighbourhood V(',num2str(ic),',',num2str(jc),',',num2str(kc),'),',tname,'=', num2str(val) ) );

% print only labels of layers where 'val' was found
iv = num2str( sort(unique( ID(:,1) ) )' );
jv = num2str( sort(unique( ID(:,2) ) )' );
kv = num2str( sort(unique( ID(:,3) ) )' ); 

xlabel( strcat('$ J=[',jv,'] $'),'interpreter','latex' );
ylabel( strcat('$ I=[',iv,'] $'),'interpreter','latex' );
zlabel( strcat('$ K=[',kv,'] $'),'interpreter','latex' );

hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',alpha);
axis equal; view(3); axis tight; axis vis3d; grid off;
set(gca,'XTick',[],'YTick',[],'ZTick',[])
colormap(cMap); 
caxis( [ min(PV(:)), max(PV(:)) ] );
%colorbar

if p
    
    % print to file    
    print('-depsc2','-r0',fullfile( '../figs/', ...
    strcat('VoxelNeighByValue_IC_',num2str( ic ),...
                            '_JC_',num2str( jc ),...
                            '_KC_',num2str( kc ),...
                            '_P_' ,num2str( P  ),...
                            '_',tname,'_'       ,...
                                   num2str(val ) ) ) );  
end


end

