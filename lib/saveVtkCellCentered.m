function saveVtkCellCentered( M, fname, cname, sname )
% SAVEVTKCELLCENTERED export 3D data array to VTK 
%                     in voxel mode.
%
%       input:  
%               M: 3D data array
%           fname: file name 
%           cname: cell_data name
%           sname: scalars name
%
%       Remark: size(M) is considered the number of voxels. 
%               Therefore, Paraview will represent a structured
%               grid with 'nverts' number of vertices
%
%   -------------
%   Dr. Gustavo Peixoto, 
%   Acks to #utkarshayachit

% checking
if ndims(M) ~= 3
    error('Input argument is not a 3D array.');
end

% grid size
I = size(M,1);
J = size(M,2);
K = size(M,3);

nvxs = I*J*K;                % number of voxels
nverts =  (I+1)*(J+1)*(K+1); % number of vertices

% file path
f = strcat('../vtk/',fname,'.vtk');
fid = fopen(f,'w');

% header 
fprintf(fid,['# vtk DataFile Version 4.0\n' ...
            'RESERVOIR CHARACTERIZATION \n']); 
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET STRUCTURED_POINTS');
fprintf(fid,'%s\n',strcat(['DIMENSIONS',' ',num2str(I+1,'%u'),' ',num2str(J+1,'%u'),' ', num2str(K+1,'%u') ]));
fprintf(fid,'%s\n',strcat('SPACING 1 1 1'));
fprintf(fid,'%s\n',strcat('ORIGIN 1 1 1'));
fprintf(fid,'%s %u \n','CELL_DATA',nvxs);
fprintf(fid,'%s\n',strcat(['SCALARS',' ',cname,' ','float']));
fprintf(fid,'%s\n',strcat('LOOKUP_TABLE default'));

% writing voxel data
v = reshape( M, [ numel(M) 1 ] );
for i = 1:length(v)
    fprintf(fid,'%f\n',v(i));
end

% vertice data
fprintf(fid,'%s %u \n','POINT_DATA',nverts);
fprintf(fid,'%s\n',strcat(['SCALARS',' ',sname,' ','float']));
fprintf(fid,'%s\n','LOOKUP_TABLE default');

%-------------- "Dummy" point data 
% Strategy: create a new matrix by extending the original one
%           by 1 slice, by duplicating the last one.

Maux = zeros( I+1,J+1,K+1 );    % auxiliary matrix
for k = 1:K    
    Mk = M(1:end,1:end,k);
    Mki = Mk(end,1:end);        % last I - row
    Mkj = Mk(1:end,end);        % last J - column
    Mkj = [ Mkj; Mkj(end) ];
    aux = cat(1,Mk,Mki);        % concatenate I
    aux = cat(2,aux,Mkj);       % concatenate J
    Maux(1:end,1:end,k) = aux;    
    
end
Maux(1:end,1:end,end) = Maux(1:end,1:end,end-1); % duplicate K

% writing point data 
%  
% row 1 (k = 1), run i,j
% row 2 (k = 2), run i,j
%  .
%  .
%  .
% row K (k = K), run i,j
%
for k = 1:size(Maux,3);
    aux = Maux(1:end,1:end,k);
    aux = reshape( aux, [ 1 numel(aux) ] );    
    fprintf(fid, strcat( repmat('%f ',[ 1 length(aux) ]),'\n'), aux );
end

fclose(fid);

end

