%% mainComponent2Image.m

clear all; close all;
%% 

% load DRT matrix
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;
DRT = 0*DRT;

ic = 45; jc = 68; % well
P = [14,14]; % VOI rings
drt = 13;

aux = load('../mat/DRT_VOI_13_Well_I45_J68.mat');
VOI = aux.VOISt;

% saving dir
svdir = '../img/';

% image format 
fmt = '.jpg';

j = 2;
ind = VOI.compVoxelInds{j};
DRT(ind) = 1;
vol = DRT( ic-P(1):ic+P(1), jc-P(2):jc+P(2), : );
vol = mat2gray(vol);


for k = 1:size(vol,3);
    aux = num2str(k); 
    reldir = strcat('subvol_drt',num2str(drt),'_comp',num2str(j),'_img');
    [~,~,~] = mkdir(svdir,reldir);    
    name = strcat('subvol_drt',num2str(drt),'_comp',num2str(j),'_');    
    img = vol(:,:,k);
    imwrite(img, strcat( fullfile(svdir, reldir, strcat(name,aux) ),fmt) );        
end
    
