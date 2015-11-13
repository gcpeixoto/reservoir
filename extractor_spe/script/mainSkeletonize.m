%% mainSkeletonize.m 
%   authors: Dr. Gustavo Peixoto de Oliveira
%            Dr. Waldir Leite Roque
%            @Federal University of Paraiba
%   mail: gustavo.oliveira@ci.ufpb.br    
%   date: Nov 11th, 2015      
%
%   Description: computes the medial axis of a given DRT region
%                and delivers the coordinates of the points so as
%                to feature the skeleton of the 3D volume.
%


clear all; close all; clc;


aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;
DRT = 0*DRT;
DRT = logical(DRT);

VOI = load('../mat/DRT_VOI_8_Well_I45_J68.mat');
VOI = VOI.VOISt;
idv = VOI.compVoxelInds{1};

DRT(idv) = true;
save('MAT.mat','DRT');
DRTS = smooth3(DRT,'gaussian');
DRTS2 = smooth3(DRT,'box');

skel = Skeleton3D(DRT);
figure();
%col = [.7 .7 .8];
col = rgb('LightBlue');
hiso = patch(isosurface(DRTS,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(DRTS,0),'FaceColor',col,'EdgeColor','none');
%axis equal;axis off;
lighting phong;
isonormals(DRTS,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'s','Markersize',4,'MarkerFaceColor','k','Color','k');            
set(gcf,'Color','white');
view(140,80)

