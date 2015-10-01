function [phi,per] = loadFiles(phiname,pername)

% delete old files
disp('Deleting old files in the tree...');
delete('extractorSPE.log', ...
       '../figs/*.*', ...       
       '../vtk/*.*' , ...
       '../csv/*.*' , ...
       '../csv/condensed/*.*' , ...
       '../mat/*.*');

disp('Loading files...');

fid = fopen(phiname);
fjd = fopen(pername);
if fopen(fid) == -1
    error('Porosity file not found.');
end
if fopen(fjd) == -1
    error('Permeability file not found.');
end

phi = load(phiname,'-ascii');
per = load(pername,'-ascii');
disp('Files loaded.');   

end

