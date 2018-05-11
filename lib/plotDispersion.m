function plotDispersion( ws,wMat,scname,ia,ja,i,mt,mfc,varargin)
%PLOTDISPERSION scatter plot of dispersion of the scalar along the well.

% input:
%          ws: well scalar field 
%        wMat: well z-coords
%      scname: scalar name ('phi','kx','ky','kz');
%          ia: well surface I-coord
%          ja: well surface J-coord
%           i: loop index
%          mt: marker type (default: black)
%         mfc: marker face color (default: black)
%       
% Use: 
%   plotDispersion( wphi{i},wMat{i,3},'phi',ia,ja,i,'d','r');
%
%
% 
% Remark: fliplr and 'Ydir reverse' are necessary to view 
%         the dispersion data of the reservoir's with z = 0 at the
%         free surface and z = zmax at the bottom.        
% 
%   depth
%   z @surface (voxel 1)
%     | o
%     |     o
%     |    
%     |o   
%     |  o 
%     |       o
%      ---------> phi, k(.)
%          

if nargin == 6
    mt = 'k'; mfc = mt; 
end

% graph padding
dpx = 0.3; 
dpy = 2;        
    
scatter(ws,fliplr(wMat),'fill',mt,'MarkerFaceColor',mfc)
h = gca;
set(h,'YGrid','on','YDir','reverse');
xlim( [ (1.0 - dpx)*min(ws) (1.0 + dpx)*max(ws) ] )
ylim( [ -dpy length(wMat)+dpy ] );        
title( strcat('Well (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
ylabel('z'); 

switch scname
    case 'phi'
        xlabel('$ \phi $','interpreter','latex');
        
    case 'kx'
        xlabel('$ \kappa_x $','interpreter','latex');
        
    case 'ky'
        xlabel('$ \kappa_y $','interpreter','latex');
        
    case 'kz'
        xlabel('$ \kappa_z $','interpreter','latex');
end


end

