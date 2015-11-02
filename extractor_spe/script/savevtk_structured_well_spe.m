function savevtk_structured_well_spe( count, I, J, phi, kx, ky, kz, name )
%  SAVEVTK_STRUCTURED_WELL_SPE Plot each well as a column of voxels.
%     input: 
%     count: well counter
%         I: well surface coordinate
%         J: well surface coordinate
%       phi: 3D array for porosity 
%      k(.): 3D array for permeabilities 
%      name: file name
%     
%     output: 
%         .vtk file
   

i = 2; % fixed voxel sizes
j = 2;
if ( length(phi) ~= length(kx) || ... 
     length(kx) ~= length(ky)  || ... 
     length(ky) ~= length(kz)  || ...
     length(kz) ~= length(phi)  )
     error('Vectors must have the same length.');
else
    k = length(phi);
end

fprintf('Processing VTK for random well %d (%d,%d)... \n',count,I,J);

fname = strcat(name,'.vtk');
fid = fopen(fname,'w');

fprintf(fid,['# vtk DataFile Version 2.0\n' ...
            'spe_data\n']); 
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET STRUCTURED_POINTS');
fprintf(fid,'%s\n',strcat(['DIMENSIONS',' ',num2str(i,'%u'),' ',num2str(j,'%u'),' ', num2str(k,'%u') ]));
fprintf(fid,'%s\n','ASPECT_RATIO 1 1 1');
fprintf(fid,'%s\n','ORIGIN 1 1 1');
fprintf(fid,'%s\n',strcat(['POINT_DATA',' ',num2str(i*j*k,'%u')]));

fprintf(fid,'%s\n',strcat('SCALARS phi_',num2str(I),'_',num2str(J),' float 1'));
fprintf(fid,'%s\n','LOOKUP_TABLE default');
% vector entries
for i = 1:length(phi)
    for j = 1:4
        fprintf(fid,'%g \n',phi(i));
    end
end

fprintf(fid,'%s\n',strcat('SCALARS kx_',num2str(I),'_',num2str(J),' float 1'));
fprintf(fid,'%s\n','LOOKUP_TABLE default');
% vector entries
for i = 1:length(kx)
    for j = 1:4
        fprintf(fid,'%g \n',kx(i));
    end
end

fprintf(fid,'%s\n',strcat('SCALARS ky_',num2str(I),'_',num2str(J),' float 1'));
fprintf(fid,'%s\n','LOOKUP_TABLE default');
% vector entries
for i = 1:length(ky)
    for j = 1:4
        fprintf(fid,'%g \n',ky(i));
    end
end
fprintf(fid,'%s\n',strcat('SCALARS kz_',num2str(I),'_',num2str(J),' float 1'));
fprintf(fid,'%s\n','LOOKUP_TABLE default');
% vector entries
for i = 1:length(kz)
    for j = 1:4
        fprintf(fid,'%g \n',kz(i));
    end
end

disp('VTK files saved.')

fclose(fid);

end

