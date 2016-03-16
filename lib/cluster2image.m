function cluster2image( DRT,drtVal,ilims,jlims,klims,fname )
%CLSUTER2IMAGE export cluster data to image sequence

% cuboid grid formed by all voxels in the cluster's bounding box
I = ilims(1):ilims(2);
J = jlims(1):jlims(2);
K = klims(1):klims(2);

DRT = DRT(I,J,K);

DRT = DRT == drtVal;
DRTclust = mat2gray(DRT);

% saving dir
svdir = '../img/';
    [~,~,~] = mkdir(svdir,fname); % creates dir       
    
for k = 1:length(K)
    img = DRTclust(:,:,k);
    imwrite(img, strcat( fullfile(svdir, fname,strcat('z',num2str(K(k)))),'.jpg') );        
end


