function plotMetricField( cvc,mfield,mname,idComp,val,cvcfar,sz,fmt,varargin )
% PLOTMETRICFIELD plots the distribution of centrality measures 
%                 along the connected component
% input:
%          cvc: cell with the voxel coordinates for a specified component
%       mfield: cell with the centrality values
%        mname: centrality metrics ('deg','bet', or 'clo')
%       idComp: connected component id
%          val: DRT value
%       cvcfar: coordinates of specific points (see mainVOIMetricsAnalyzer)
%           sz: plot marker size
%          fmt: format for figure printing 

if nargin == 5
    cvcfar = cvc;
    sz = 40;
    sz2 = sz;
    fmt = 'pdf';        
elseif nargin == 6
    sz = 40;
    sz2 = 300;
    fmt = 'pdf';        
end

figure    
scatter3( cvc(:,2), cvc(:,1), cvc(:,3),sz,mfield,'filled');
hold on
scatter3(cvcfar(:,2),cvcfar(:,1),cvcfar(:,3),sz2,'k'); % highlight farthest points 
grid off
axis equal; view(3); axis tight; axis vis3d; grid off;             
set(gca,'ZDir','reverse');    

if strcmp(mname,'deg')      % degree centrality
    flag = '\delta';
elseif strcmp(mname,'bet')  % betweeness centrality
    flag = '\beta';
elseif strcmp(mname,'clo')  % closeness centrality
    flag = '\gamma';
end    
straux = strcat('$; \, ',flag,' $');
str = strcat( '$ C $',num2str(idComp),'$ - DRT=$ ', num2str( val ),straux);
title( str ,'interpreter','latex');
xlabel('J');
ylabel('I');
zlabel('K');    
colorbar

% print
if strcmp(fmt,'eps') 
    print('-depsc2','-r0',strcat( '../figs/graphPath/', ...
        'Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ),'_',mname,'.eps' ) );  
    elseif strcmp(fmt,'pdf') 
    print('-dpdf','-r0',strcat( '../figs/graphPath/', ...
        'Voxel_Graph_Component',num2str(idComp),'_DRT_',num2str( val ),'_',mname,'.pdf' ) );  
end  



