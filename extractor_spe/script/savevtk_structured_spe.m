function savevtk_structured_spe( i, j, k, phi, kx, ky, kz, name )
%  SAVEVTK_STRUCTURED_SPE Plot structured mesh of the full reservoir.
%     input:     
%         i,j,k: reservoir sizes (length, width, depth)
%       phi: 3D array for porosity 
%      k(.): 3D array for permeabilities 
%      name: file name
%     
%     output: 
%         .vtk file
   

if ( size(phi,1) > 1 && size(phi,2) > 1 && size(phi,3) > 1 && ...
     size(kx,1)  > 1 && size(kx,2)  > 1 && size(kx,3)  > 1 && ... 
     size(ky,1)  > 1 && size(ky,2)  > 1 && size(ky,3)  > 1 && ...
     size(kz,1)  > 1 && size(kz,2)  > 1 && size(kz,3)  > 1  )
   
        disp('Converting 3D matrices to vector for VTK plotting...');
        phi = reshape( phi, [ numel(phi) 1 ] );
        kx = reshape( kx, [ numel(kx) 1 ] );
        ky = reshape( ky, [ numel(ky) 1 ] );
        kz = reshape( kz, [ numel(kz) 1 ] );

end

fname = strcat(name,'.vtk');
fid = fopen(fname,'w');

fprintf(fid,['# vtk DataFile Version 2.0\n' ...
            'spe_data\n']); 
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET STRUCTURED_POINTS');
fprintf(fid,'%s\n',strcat(['DIMENSIONS',' ',num2str(i,'%u'),' ',num2str(j,'%u'),' ', num2str(k,'%u') ]));
fprintf(fid,'%s\n',strcat('ASPECT_RATIO 1 1 1'));
fprintf(fid,'%s\n',strcat('ORIGIN 1 1 1'));
fprintf(fid,'%s\n',strcat(['POINT_DATA',' ',num2str(i*j*k,'%u')]));

fprintf(fid,'%s\n',strcat('SCALARS porosity float 1'));
fprintf(fid,'%s\n',strcat('LOOKUP_TABLE default'));
% vector entries
for i = 1:length(phi)
    fprintf(fid,'%g \n',phi(i));
end

fprintf(fid,'%s\n',strcat('SCALARS kx float 1'));
fprintf(fid,'%s\n',strcat('LOOKUP_TABLE default'));
% vector entries
for i = 1:length(kx)
    fprintf(fid,'%g \n',kx(i));
end

fprintf(fid,'%s\n',strcat('SCALARS ky float 1'));
fprintf(fid,'%s\n',strcat('LOOKUP_TABLE default'));
% vector entries
for i = 1:length(ky)
    fprintf(fid,'%g \n',ky(i));
end

fprintf(fid,'%s\n',strcat('SCALARS kz float 1'));
fprintf(fid,'%s\n',strcat('LOOKUP_TABLE default'));
% vector entries
for i = 1:length(kz)
    fprintf(fid,'%g \n',kz(i));
end


fclose(fid);

end

