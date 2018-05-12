%% mainProcessClusterFit
% Determines locations to setup nonuniform 5-spot well patterns based 
% on the best-fit ellipsoid obtained from the convex hull points 
% of the selected cluster
%
% METHODOLOGY
% ===========
%
% - Consider G a cluster with arbitrary shape and c the 
%   voxels (vertices) that belong to the convex hull H of G 
%
% - The best-fit ellipsoid based on c points returns 
%   the 3 main direction vectors, the 3 radii a,b,c of each axis 
%   and the center point 
%
%      o X1                      X1: fit-point (+x)                
%                                X2: fit-point (-x)
%        P1-c                    Y1: fit-point (+y)                
%        |  |       o Y1         Y2: fit-point (-y)
%        c  c--c-P2              Z1: fit-point (+z) (not used for vertical)
%        |    C   |              Z2: fit-point (-z) (not used for vertical)
%  o Y2  P4-c  c--c
%           |  |                  C: center
%           c-P3
%                                Pk are such that min {d(Xk,H)} in L2-norm
%                 o X2
%  RELATIONS
%  ========
%
%  X1 = C + a*E1        X2 = C - a*E1
%  Y1 = C + b*E2        Y2 = C - b*E2
%  Z1 = C + c*E3        Z2 = C - c*E3, for orthonormal basis {E1,E2,E3}
%
%  NONUNIFORM 5-SPOT PATTERNS
%  ==========================
% 
%  - Formed by the PRODUCER column placed at the max closeness 
%    vertex M = (Mx,My,Mz) z-span, i.e. (Mx,My,k), k \in klims.
%    + 4 INJECTOR columns placed at 
%    Pj = (Pj_x,Pj_y,k), k \in klims, j = 1,2,3,4. 
%  
%  ROTATED CONFIGURATIONS
%  ======================
%
%  - This script also computes this 5-spot pattern rotatated 
%    by an angle theta in relation to the reference frame 
%    {C,X1<->X2,Y1<->Y2}
%% Load cluster info

clear, close all;  
load('../mat/Field/DRT_13.mat','drtSt');
load('../mat/Field/DRT_13_MetricsData.mat','metrics')
load('../mat/Field/DRT_13_LinRegrData.mat','linregr')

c = 2; % cluster ID
theta = pi/3; % angle (in radians)

% voxels
cvc = drtSt.compVoxelCoords{c}; 

% cluster limits
iG = cvc(:,1); jG = cvc(:,2); kG = cvc(:,3);

% cluster bounding box limits
im=min(cvc(:,1)); iM=max(cvc(:,1)); ilims = [im,iM];
jm=min(cvc(:,2)); jM=max(cvc(:,2)); jlims = [jm,jM];
km=min(cvc(:,3)); kM=max(cvc(:,3)); klims = [km,kM];

% max closeness value, coords and index
maxc = max(metrics.closenessCentrality{c}); 
iMaxc = find(metrics.closenessCentrality{c} == maxc);
vcMaxc = cvc(iMaxc,:);
MCX = vcMaxc(1); MCY = vcMaxc(2); MCZ = vcMaxc(3);

% column neighbours inside the cluster z-range (maxC)
zz = unique(kG); 
lz = length(zz);
cols = [ ones(lz,1)*vcMaxc(1) ones(lz,1)*vcMaxc(2) zz ]; 

% Find convex hull 
k = convhull(cvc);
ind = unique(k(:));
vch = [iG(ind),jG(ind),kG(ind)]; % coords

% ellipsoid fit over convex hull points
[C,R,E,P] = ellipsoid_fit(vch);

% ellipsoid extrema points over principal axes
P1 = C + R(1)*E(:,1); % +X
P2 = C + R(2)*E(:,2); % +Y
P3 = C + R(3)*E(:,3); % +Z

P4 = C - R(1)*E(:,1); % -X
P5 = C - R(2)*E(:,2); % -Y
P6 = C - R(3)*E(:,3); % -Z

% store
clusterFitData.drtValue = drtSt.value;
clusterFitData.ncomp = c;
clusterFitData.iMaxc = iMaxc;
clusterFitData.vcMaxc = vcMaxc;
clusterFitData.colNeighsMaxC= cols;
clusterFitData.iHull = ind;
clusterFitData.vcHull = vch;
clusterFitData.ellipsoidFitCenter = C;
clusterFitData.ellipsoidFitRadii = R;
clusterFitData.ellipsoidFitEvecs = E;
clusterFitData.ellipsoidFitParams = P;
clusterFitData.rotationAngleDeg = 180*theta/pi;

% vcMaxc - local coords
aux = global2LocalCluster(vcMaxc(1),vcMaxc(2),vcMaxc(3),im,jm,km);
clusterFitData.vcMaxcLocalCoords = aux;

%% Find 6 farthest points in the cluster from ellipsoid fit points

% directions
pos = {'+X','-X','+Y','-Y','+Z','-Z'};

% ----------
% Rotation matrix in relation to an axis
% - matrix:
%      R = cos(theta)I + sin(theta)UX + (1-cos(theta))*[u dyad u]
% - angle:
%      theta
% - axis direction unit vector: 
%     (ux,uy,uz)
% ----------

ux = E(1,3); uy = E(2,3); uz = E(3,3); % direction vector

% matrices
UX = [ 0 -uz uy; 
       uz 0 -ux;
      -uy ux 0 ];
  
UU = [ ux*ux ux*uy ux*uz; 
       ux*uy uy*uy uy*uz;
       ux*uz uy*uz uz*uz ];

% rotation matrix
ROT = cos(theta)*eye(3) + sin(theta)*UX + (1-cos(theta))*UU;

% xy
[RC, RP1, RP2, RP4, RP5] = deal(ROT*C, ROT*P1, ROT*P2, ROT*P4, ROT*P5);

% displacing points in relation to retake center
D = RC - C;
RP1 = RP1 - D;
RP2 = RP2 - D;
RP4 = RP4 - D;
RP5 = RP5 - D;

% distances 
dc = zeros(size(vch,1),1);
dcr = zeros(size(vch,1),1); % rotated 
for p = 1:length(pos)

    posi = pos{p};
    
    switch posi
        
        case '+X'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P1 - P' ).^2 ) );                                        
                dcr(i) = sqrt( sum( ( RP1 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);
            
            [dcr,auxr] = sortrows(dcr);
            ipr = auxr(1);
            vcr = vch(ipr,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 
            
            %---- SAVE
            
            % save colNeighsX1
            clusterFitData.colNeighsX1 = col_neighbours;
            
            % save local coords - colNeighsX1 
            aux = global2LocalCluster(col_neighbours(:,1),...
                                      col_neighbours(:,2),...
                                      col_neighbours(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsX1LocalCoords = aux;
            
            % save index X1
            clusterFitData.ipX1 = ip;                        
            
            % save voxel X1
            clusterFitData.vcX1 = vc;
            
            % save local coords - vcX1 
            aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
            clusterFitData.vcX1LocalCoords = aux;
            
            % save colNeighsX1R
            clusterFitData.colNeighsX1R = col_neighboursR;
            
            % local coords - colNeighsX1R 
            aux = global2LocalCluster(col_neighboursR(:,1),...
                                      col_neighboursR(:,2),...
                                      col_neighboursR(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsX1RLocalCoords = aux;
            
            % index X1R
            clusterFitData.ipX1R = ipr;
            
            % voxel X1R
            clusterFitData.vcX1R = vcr;
            
            % local coords - vcR
            aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
            clusterFitData.vcX1RLocalCoords = aux;
            

        case '-X'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P4 - P' ).^2 ) );                                        
                dcr(i) = sqrt( sum( ( RP4 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);
            
            [dcr,auxr] = sortrows(dcr);
            ipr = auxr(1);
            vcr = vch(ipr,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 
            
            %---- SAVE
            
            % save colNeighsX2
            clusterFitData.colNeighsX2 = col_neighbours;
            
            % save local coords - colNeighsX2
            aux = global2LocalCluster(col_neighbours(:,1),...
                                      col_neighbours(:,2),...
                                      col_neighbours(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsX2LocalCoords = aux;
            
            % save index X2
            clusterFitData.ipX2 = ip;                        
            
            % save voxel X2
            clusterFitData.vcX2 = vc;
            
            % save local coords - vcX2
            aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
            clusterFitData.vcX2LocalCoords = aux;
            
            % save colNeighsX2R
            clusterFitData.colNeighsX2R = col_neighboursR;
            
            % local coords - colNeighsX2R 
            aux = global2LocalCluster(col_neighboursR(:,1),...
                                      col_neighboursR(:,2),...
                                      col_neighboursR(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsX2RLocalCoords = aux;
            
            % index X2R
            clusterFitData.ipX2R = ipr;
            
            % voxel X2R
            clusterFitData.vcX2R = vcr;
            
            % local coords - vcR
            aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
            clusterFitData.vcX2RLocalCoords = aux;
      
        case '+Y'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P2 - P' ).^2 ) );     
                dcr(i) = sqrt( sum( ( RP2 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);
            
            [dcr,auxr] = sortrows(dcr);
            ipr = auxr(1);
            vcr = vch(ipr,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 
            
            %---- SAVE
            
            % save colNeighsY1
            clusterFitData.colNeighsY1 = col_neighbours;
            
            % save local coords - colNeighsY1 
            aux = global2LocalCluster(col_neighbours(:,1),...
                                      col_neighbours(:,2),...
                                      col_neighbours(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsY1LocalCoords = aux;
            
            % save index Y1
            clusterFitData.ipY1 = ip;                        
            
            % save voxel Y1
            clusterFitData.vcY1 = vc;
            
            % save local coords - vcY1 
            aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
            clusterFitData.vcY1LocalCoords = aux;
            
            % save colNeighsY1R
            clusterFitData.colNeighsY1R = col_neighboursR;
            
            % local coords - colNeighsY1R 
            aux = global2LocalCluster(col_neighboursR(:,1),...
                                      col_neighboursR(:,2),...
                                      col_neighboursR(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsY1RLocalCoords = aux;
            
            % index Y1R
            clusterFitData.ipY1R = ipr;
            
            % voxel Y1R
            clusterFitData.vcY1R = vcr;
            
            % local coords - vcR
            aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
            clusterFitData.vcY1RLocalCoords = aux;
    
        case '-Y'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P5 - P' ).^2 ) );  
                dcr(i) = sqrt( sum( ( RP5 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);
            
            [dcr,auxr] = sortrows(dcr);
            ipr = auxr(1);
            vcr = vch(ipr,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            col_neighboursR = [ ones(lz,1)*vcr(1) ones(lz,1)*vcr(2) zz ]; 
            
            %---- SAVE
            
            % save colNeighsY2
            clusterFitData.colNeighsY2 = col_neighbours;
            
            % save local coords - colNeighsY2
            aux = global2LocalCluster(col_neighbours(:,1),...
                                      col_neighbours(:,2),...
                                      col_neighbours(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsY2LocalCoords = aux;
            
            % save index Y2
            clusterFitData.ipY2 = ip;                        
            
            % save voxel Y2
            clusterFitData.vcY2 = vc;
            
            % save local coords - vcY2
            aux = global2LocalCluster(vc(1),vc(2),vc(3),im,jm,km);
            clusterFitData.vcY2LocalCoords = aux;
            
            % save colNeighsY2R
            clusterFitData.colNeighsY2R = col_neighboursR;
            
            % local coords - colNeighsY2R 
            aux = global2LocalCluster(col_neighboursR(:,1),...
                                      col_neighboursR(:,2),...
                                      col_neighboursR(:,3),...
                                      im,jm,km);
            clusterFitData.colNeighsY2RLocalCoords = aux;
            
            % index Y2R
            clusterFitData.ipY2R = ipr;
            
            % voxel Y2R
            clusterFitData.vcY2R = vcr;
            
            % local coords - vcR
            aux = global2LocalCluster(vcr(1),vcr(2),vcr(3),im,jm,km);
            clusterFitData.vcY2RLocalCoords = aux;
            
        case '+Z'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P3 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            
            clusterFitData.colNeighsZ1 = col_neighbours;
            clusterFitData.ipZ1 = ip;
            clusterFitData.vcZ1 = vc;
            
        case '-Z'
            
            % span convex hull
            for i = 1:size(vch,1) 
                P = vch(i,:);       
                dc(i) = sqrt( sum( ( P6 - P' ).^2 ) );                                        
            end

            % sort distances and get minimum one
            [dc,aux] = sortrows(dc);
            ip = aux(1);
            vc = vch(ip,:);

            % column neighbours inside the cluster z-range
            zz = unique(kG); 
            lz = length(zz);
            col_neighbours = [ ones(lz,1)*vc(1) ones(lz,1)*vc(2) zz ]; 
            
            clusterFitData.colNeighsZ2 = col_neighbours;
            clusterFitData.ipZ2 = ip;
            clusterFitData.vcZ2 = vc;
            
    end
end

% cluster normal 5-spot perforations
zz = unique(kG); 
lz = length(zz);
col_neighbours_INJ1 = [ ones(lz,1)*min(iG) ones(lz,1)*min(jG) zz ]; 
col_neighbours_INJ2 = [ ones(lz,1)*max(iG) ones(lz,1)*min(jG) zz ]; 
col_neighbours_INJ3 = [ ones(lz,1)*max(iG) ones(lz,1)*max(jG) zz ]; 
col_neighbours_INJ4 = [ ones(lz,1)*min(iG) ones(lz,1)*max(jG) zz ]; 
PM  = round(mean(cvc));
col_neighbours_PROD = [ ones(lz,1)*PM(1) ones(lz,1)*PM(2) zz ]; 

 
%% Visualization 

% cluster points 
%scatter3(iG,jG,kG,'k'),'filled');
axis equal vis3d on
hold on

% ==== ORIGINAL FRAME

% original center point 
scatter3(C(1),C(2),C(3),'m','filled');

% extrema points +
scatter3(P1(1),P1(2),P1(3),'r','filled');
scatter3(P2(1),P2(2),P2(3),'g','filled');
scatter3(P3(1),P3(2),P3(3),'b','filled');

% extrema points -
scatter3(P4(1),P4(2),P4(3),'r','filled');
scatter3(P5(1),P5(2),P5(3),'g','filled');
scatter3(P6(1),P6(2),P6(3),'b','filled');

% rotated points (XY direction)
scatter3(RP1(1),RP1(2),RP1(3),'dr','filled');
scatter3(RP2(1),RP2(2),RP2(3),'dg','filled');

scatter3(RP4(1),RP4(2),RP4(3),'dr','filled');
scatter3(RP5(1),RP5(2),RP5(3),'dg','filled');

% max closeness point
%scatter3(vcMaxc(1),vcMaxc(2),vcMaxc(3),'y','filled');

% convex hull
%scatter3(vch(:,1),vch(:,2),vch(:,3),'c','filled');

% principal axes (original ellipsoid)
line([P1(1),P4(1)],[P1(2),P4(2)],[P1(3),P4(3)])
line([P2(1),P5(1)],[P2(2),P5(2)],[P2(3),P5(3)])
line([P3(1),P6(1)],[P3(2),P6(2)],[P3(3),P6(3)])

% rotated axes
line([RP1(1),RP4(1)],[RP1(2),RP4(2)],[RP1(3),RP4(3)])
line([RP2(1),RP5(1)],[RP2(2),RP5(2)],[RP2(3),RP5(3)])

% chosen points for injection + associated column 

% maxC column
scatter3(clusterFitData.colNeighsMaxC(:,1),... 
         clusterFitData.colNeighsMaxC(:,2),... 
         clusterFitData.colNeighsMaxC(:,3),... 
         'c','filled');

% adjacent columns      
scatter3(clusterFitData.colNeighsX1(:,1),... 
         clusterFitData.colNeighsX1(:,2),... 
         clusterFitData.colNeighsX1(:,3),... 
         'm','filled');

scatter3(clusterFitData.colNeighsX2(:,1),... 
         clusterFitData.colNeighsX2(:,2),... 
         clusterFitData.colNeighsX2(:,3),... 
         'm','filled');     
     
scatter3(clusterFitData.colNeighsY1(:,1),... 
         clusterFitData.colNeighsY1(:,2),... 
         clusterFitData.colNeighsY1(:,3),... 
         'm','filled');


scatter3(clusterFitData.colNeighsY2(:,1),... 
         clusterFitData.colNeighsY2(:,2),... 
         clusterFitData.colNeighsY2(:,3),... 
         'm','filled');     

scatter3(clusterFitData.colNeighsZ1(:,1),... 
         clusterFitData.colNeighsZ1(:,2),... 
         clusterFitData.colNeighsZ1(:,3),... 
         'm','filled');


scatter3(clusterFitData.colNeighsZ2(:,1),... 
         clusterFitData.colNeighsZ2(:,2),... 
         clusterFitData.colNeighsZ2(:,3),... 
         'm','filled');     
        
% adjacent rotated columns      
scatter3(clusterFitData.colNeighsX1R(:,1),... 
         clusterFitData.colNeighsX1R(:,2),... 
         clusterFitData.colNeighsX1R(:,3),... 
         'dm');

scatter3(clusterFitData.colNeighsX2R(:,1),... 
         clusterFitData.colNeighsX2R(:,2),... 
         clusterFitData.colNeighsX2R(:,3),... 
         'dm');     
     
scatter3(clusterFitData.colNeighsY1R(:,1),... 
         clusterFitData.colNeighsY1R(:,2),... 
         clusterFitData.colNeighsY1R(:,3),... 
         'dm');

scatter3(clusterFitData.colNeighsY2R(:,1),... 
         clusterFitData.colNeighsY2R(:,2),... 
         clusterFitData.colNeighsY2R(:,3),... 
         'dm');            
             
% points
scatter3(clusterFitData.vcX1(1),... 
         clusterFitData.vcX1(2),... 
         clusterFitData.vcX1(3),... 
         'y','filled');

scatter3(clusterFitData.vcX2(1),... 
         clusterFitData.vcX2(2),... 
         clusterFitData.vcX2(3),... 
         'y','filled');

scatter3(clusterFitData.vcY1(1),... 
         clusterFitData.vcY1(2),... 
         clusterFitData.vcY1(3),... 
         'y','filled');

scatter3(clusterFitData.vcY2(1),... 
         clusterFitData.vcY2(2),... 
         clusterFitData.vcY2(3),... 
         'y','filled');
     
scatter3(clusterFitData.vcZ1(1),... 
         clusterFitData.vcZ1(2),... 
         clusterFitData.vcZ1(3),... 
         'y','filled');

scatter3(clusterFitData.vcZ2(1),... 
         clusterFitData.vcZ2(2),... 
         clusterFitData.vcZ2(3),... 
         'y','filled');    

% rotated points
scatter3(clusterFitData.vcX1R(1),... 
         clusterFitData.vcX1R(2),... 
         clusterFitData.vcX1R(3),... 
         'dy');

scatter3(clusterFitData.vcX2R(1),... 
         clusterFitData.vcX2R(2),... 
         clusterFitData.vcX2R(3),... 
         'dy');

scatter3(clusterFitData.vcY1R(1),... 
         clusterFitData.vcY1R(2),... 
         clusterFitData.vcY1R(3),... 
         'dy');

scatter3(clusterFitData.vcY2R(1),... 
         clusterFitData.vcY2R(2),... 
         clusterFitData.vcY2R(3),... 
         'dy');
     
              
% ==== DISPLACED FRAME
%{
% displacement line 
line([C(1),MCX],[C(2),MCY],[C(3),MCZ])


% displaced frame 
dx = C(1) - MCX;
dy = C(2) - MCY;
dz = C(3) - MCZ;

D = [dx,dy,dz]';

P7 = P1 - D;
P8 = P2 - D;
P9 = P3 - D;

P10 = P4 - D;
P11 = P5 - D;
P12 = P6 - D;

% extrema points - displaced
scatter3(P7(1),P7(2),P7(3),'r','filled', 'marker','s');
scatter3(P8(1),P8(2),P8(3),'g','filled', 'marker','s');
scatter3(P9(1),P9(2),P9(3),'b','filled', 'marker','s');

scatter3(P10(1),P10(2),P10(3),'r','filled', 'marker','s');
scatter3(P11(1),P11(2),P11(3),'g','filled', 'marker','s');
scatter3(P12(1),P12(2),P12(3),'b','filled', 'marker','s');

% principal axes (displaced ellipsoid)
line([P7(1),P10(1)],[P7(2),P10(2)],[P7(3),P10(3)])
line([P8(1),P11(1)],[P8(2),P11(2)],[P8(3),P11(3)])
line([P9(1),P12(1)],[P9(2),P12(2)],[P9(3),P12(3)])
%} 

% ==== BOUNDARY TOPOLOGY  


% convex hull with shrink factor (boundary function)

CHN = boundary(cvc,1);
trisurf(CHN,cvc(:,1),cvc(:,2),cvc(:,3),'FaceColor',0.5*[1,1,1],...
    'EdgeColor',0.3*[1,1,1],...
    'FaceLighting', 'gouraud')
shading interp, axis equal


%% Draw Ellipsoid
% mind = min( [ iG jG kG ] );
% maxd = max( [ iG jG kG ] );
% nsteps = 50;
% step = ( maxd - mind ) / nsteps;
% [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
% 
% Ellipsoid = P(1) *x.*x +   P(2) * y.*y + P(3) * z.*z + ...
%           2*P(4) *x.*y + 2*P(5)*x.*z + 2*P(6) * y.*z + ...
%           2*P(7) *x    + 2*P(8)*y    + 2*P(9) * z;
% 
% p = patch( isosurface( x, y, z, Ellipsoid, -P(10) ) );
% surf(x,y,z,Ellipsoid)
% hold off;
% set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
% view( -70, 40 );
% axis vis3d equal;
% camlight;
% lighting phong;

if ~exist('../mat/clusterFitData','dir')
    mkdir('../mat/clusterFitData')
end
nm = strcat('../mat/clusterFitData/','DRT_',num2str(drtSt.value),'_C_',num2str(c),'_angle_',num2str(180*theta/pi),'_FitData.mat');
save(nm,'clusterFitData');
            
