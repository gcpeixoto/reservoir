function [ D2CFar, CVCFar, ilims, jlims, klims ] = getDists2Point( cvc, p, N )
% GETDISTS2POINT computes the distances and coordinates of 
%                points taken in a 3D bounding box to 
%                a given 3D point. 
%
%   input: 
%           cvc: coordinates of all the voxels (points)
%                of the cloud (cluster): (nx3)                
%             p: 3D base-voxel (point): (1x3)
%             N: norm exponent (1,2,P,Inf,-Inf): 'help norm'
% 
%   Motivation: function created to compute the distances 
%               from the farthest points to the voxel of
%               maximum closeness centrality in a cluster,
%               which are at the bounding box slices 
%              

if iscell(p)    
    p = p{1};
end

% cluster bounding box limits
im=min(cvc(:,1)); iM=max(cvc(:,1)); ilims = [im,iM];
jm=min(cvc(:,2)); jM=max(cvc(:,2)); jlims = [jm,jM];
km=min(cvc(:,3)); kM=max(cvc(:,3)); klims = [km,kM];

%------- voxels at bounding box slice
% - i  
indim =  cvc(:,1) == im ;
cvcim = cvc(indim,:);

% + i
indiM =  cvc(:,1) == iM ;
cvciM = cvc(indiM,:);

% - j
indjm =  cvc(:,2) == jm ;
cvcjm = cvc(indjm,:);

% + j
indjM =  cvc(:,2) == jM ;
cvcjM = cvc(indjM,:);

% - k
indkm =  cvc(:,3) == km ;
cvckm = cvc(indkm,:);

% + k
indkM =  cvc(:,3) == kM ;
cvckM = cvc(indkM,:);

%------- distances to point p
Nim = zeros( size(cvcim,1),1 );
NiM = zeros( size(cvciM,1),1 );
Njm = zeros( size(cvcjm,1),1 );
NjM = zeros( size(cvcjM,1),1 );
Nkm = zeros( size(cvckm,1),1 );
NkM = zeros( size(cvckM,1),1 );

% slice - i
for i = 1:length(Nim)
    Nim(i) = norm( cvcim(i,:) - p, N );
end

% slice + i
for i = 1:length(NiM)
    NiM(i) = norm( cvciM(i,:) - p, N );
end

% slice - j
for i = 1:length(Njm)
    Njm(i) = norm( cvcjm(i,:) - p, N );
end

% slice + j
for i = 1:length(NjM)
    NjM(i) = norm( cvcjM(i,:) - p, N );
end

% slice - k
for i = 1:length(Nkm)
    Nkm(i) = norm( cvckm(i,:) - p, N );
end

% slice + k
for i = 1:length(NkM)
    NkM(i) = norm( cvckM(i,:) - p, N );
end

%------- max distances and voxel coordinates
% slice - i
aux = find( Nim == max(Nim) );
d2pimfar = Nim(aux);
cvcimfar = cvcim(aux,:);

% slice + i
aux = find( NiM == max(NiM) );
d2piMfar = NiM(aux);
cvciMfar = cvciM(aux,:);

% slice - j
aux = find( Njm == max(Njm) );
d2pjmfar = Njm(aux);
cvcjmfar = cvcjm(aux,:);

% slice + j
aux = find( NjM == max(NjM) );
d2pjMfar = NjM(aux);
cvcjMfar = cvcjM(aux,:);

% slice - k
aux = find( Nkm == max(Nkm) );
d2pkmfar = Nkm(aux);
cvckmfar = cvckm(aux,:);

% slice + k
aux = find( NkM == max(NkM) );
d2pkMfar = NkM(aux);
cvckMfar = cvckM(aux,:);

%------- arrays (there may be identical point in norm-2)
% 6 distances-to-point
D2CFar = [ d2pimfar;
           d2piMfar;
           d2pjmfar;
           d2pjMfar;
           d2pkmfar;
           d2pkMfar; ];

% 6 voxel coords
CVCFar = [ cvcimfar;       
           cvciMfar;
           cvcjmfar;       
           cvcjMfar;
           cvckmfar;       
           cvckMfar; ];
       

end
