%% extractorDRTGraphPath

%% DEFAULTS
clear all; close all; clc; format long;
diary('extractorDRTGraphPath.log');
diary on

setOptions;

%% LOAD FILES

% file paths
phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

[PHI,KX,KY,KZ] = loadMatFiles(phiname,kxname,kyname,kzname);

%% COMPUTATION

KN = sqrt( KX.^2 + KY.^2 + KZ.^2 );
PHIZ = PHI./(1.0 - PHI);
RQI = 0.0314*sqrt( KN./PHI );
FZI = RQI./PHIZ;

[ VB, coords, inds ] = maskVoxels( FZI,Inf );
if ~isempty( coords ) 
    for e = 1:size(coords,1)
        FZI( coords(e,1), coords(e,2), coords(e,3) ) = 0.0;
    end
end

DRT = round( 2*log( FZI ) + 10.6 );

%% PLOTTING 

pltfield = false; % plot DRT = drt for all the field?

drt = sort( unique( DRT(:) ) );
drt = drt(2:end); % removes -Inf
 
%{
    tabulate:    
    DRT value | frequency | percentage
%}
tab = tabulate( DRT(:) );
valDRT = tab(:,1);
countDRT = tab(:,2);
percDRT = tab(:,3);


for m = 1:length(drt)

    % getting DRTs
    [ VBDRT, coordsDRT, indz ] = maskVoxels( DRT,drt(m) ); 

    
    if isempty(coordsDRT)
        disp( strcat( 'No voxels for DRT = ', num2str( drt(m) ) ) ); 
        continue 
    end
    
    %{
        PATH TRACKING
        =============
    
        - stores path array; 
        - gets the anchor voxel V0^{DRT=drt}(i0,j0,k0) found in coordsDRT;
        - seeks all the V^{DRT=drt}(i,j,k) such as d( V, V0 ) <= sqrt(2), 
          i.e, 8-connected neighbours and appends the linear index to path
          by updating the anchor voxel;
    
                 [ i0 j0 k0 ]           [ ind0 ]
                 [ i1 j1 k1 ]           [ ind1 ]
                 [     .    ]           [   .  ]
          path = [     .    ]  idpath = [   .  ]
                 [     .    ]           [   .  ]
                 [ in jn kn ]           [ indn ]
            
    %}
    path = [];
    idpath = [];
    c1 = coordsDRT(1,:);
    path = [ path; c1 ];
    idpath = [ idpath; indz(1) ];

    for e = 1:size(coordsDRT,1)

        ca = coordsDRT(e,1);
        cb = coordsDRT(e,2);
        cc = coordsDRT(e,3);
                
        dist = sqrt( ( ca - c1(1) )^2 + ...
                     ( cb - c1(2) )^2 + ...
                     ( cc - c1(3) )^2 ); 
        
        if dist <= sqrt(2) && dist > 0             
            path = [ path; [ ca cb cc ] ];
            idpath = [ idpath; indz(e) ];            
            % update anchor voxel
            c1(1) = ca;
            c1(2) = cb;
            c1(3) = cc;            
        end

    end

    % only graph and paths with more than 'nvpath' voxels
    nvpath = 2;
    if size(path,1) > nvpath
        
        %------------- GRAPH PATH
        figure     
        path1 = path(:,1);
        path2 = path(:,2);
        path3 = path(:,3);        
        plot3(path2,path1,path3,'LineWidth',1.5, 'Color',[0 .4 0.54], ...
            'Marker','o', 'MarkerSize',2.5, ...
            'MarkerFaceColor','k', 'MarkerEdgeColor','k')

        axis equal, axis vis3d, grid on, view(3)
        title( strcat( 'Path - DRT ', num2str( drt(m) ) ) );
        xlabel('J');
        ylabel('I');
        zlabel('K');
        pt1 = sort( unique( path1 ) );
        pt2 = sort( unique( path2 ) );
        pt3 = sort( unique( path3 ) );
        pt1 = num2str(pt1); pt2 = num2str(pt2); pt3 = num2str(pt3);
        pt1 = cellstr(pt1); pt2 = cellstr(pt2); pt3 = cellstr(pt3);
        
        set( gca,'XTick',min(path2):1:max(path2),'XTickLabel',pt2, ...
                 'YTick',min(path1):1:max(path1),'YTickLabel',pt1, ... 
                 'ZTick',min(path3):1:max(path3),'ZTickLabel',pt3 );
        
        print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
                strcat('Path_DRT_',num2str( drt(m) ) ) ) );   
        
            
        %------------- VOXEL PATH
        figure 
        [ ~, ~, indz ] = maskVoxels( DRT,drt(m) );
        [F,V,C]=ind2patch(idpath,DRT,'v');
        
        hold on;
        patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','k','FaceAlpha',0.8);
        axis equal; view(3); axis tight; axis vis3d; grid off;        
        title( strcat( 'Path - DRT ', num2str( drt(m) ) ) );
        xlabel('J');
        ylabel('I');
        zlabel('K');        
        
        print('-depsc2','-r0',fullfile( '../figs/graphPath', ...
                strcat('Voxel_Path_DRT_',num2str( drt(m) ) ) ) ); 
            
        if pltfield == true
            plotDRTField(DRT,val);
        end
    end
end


%close all
diary off