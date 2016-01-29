function plotPField( cvc,field,idComp,val,sz,fmt,varargin )
% PLOTMETRICFIELD plots the distribution of pressure
%                 along the connected component
% input:
%          cvc: cell with the voxel coordinates for a specified component
%       field: cell with the pressure values
%       idComp: connected component id
%          val: DRT value
%           sz: plot marker size
%          fmt: format for figure printing 

if nargin == 4    
    sz = 40;    
    fmt = 'pdf';        
end

figure    
scatter3( cvc(:,2), cvc(:,1), cvc(:,3),sz,field,'filled');
grid off
axis equal; view(3); axis tight; axis vis3d; grid off;             
set(gca,'ZDir','reverse');    


straux = strcat('$; \, p $');
str = strcat( '$ C $',num2str(idComp),'$ - DRT=$ ', num2str( val ),straux);
title( str ,'interpreter','latex');
xlabel('J');
ylabel('I');
zlabel('K');    
colorbar

% print
if strcmp(fmt,'eps') 
    print('-depsc2','-r0',strcat( '../figs/graphPath/', ...
        'Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ),'_press','.eps' ) );  
    elseif strcmp(fmt,'pdf') 
    print('-dpdf','-r0',strcat( '../figs/graphPath/', ...
        'Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ),'_press','.pdf' ) );  
end  



