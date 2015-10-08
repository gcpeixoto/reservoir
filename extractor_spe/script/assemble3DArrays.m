function [ PHI,KX,KY,KZ ] = assemble3DArrays( phi,per,I,J,K )
%{
    ASSEMBLE3DARRAYS Mounts 3D arrays for porosity and permeability

    input:
        phi: reservoir's porosity matrix in CMG output
        per: reservoir's permeability matrix in CMG output
      I,J,K: reservoir's length, width, depth (voxel-based)

    output:
        PHI, KX, KY, KZ: 3D arrays storing the reservoir's 
        scalar field data.

    - REMARK
    
    This strategy is, at first glance, ONLY applicable to 
    .dat files that follow the output format of CMG software. 
    The current version is tested for the SPE Project 2 
    data arrangement. It will might be useful for other 
    files as well, but further investigation is required. 

    What is done here: the original files are 'unrolled' 
    and rearranged into 3D arrays (I,J,K), thus allowing a 
    better accessibility of the information. In the end, the reservoir
    is a 3D structure with { PHI, KX, KY, KZ }(I,J,K) evaluated 
    in the fourth dimension. 

    RESERVOIR SCHEME
    ================

          K
         /
        /____________
       /|            |
      /_|__________  |
     /|            | |
    /_|__________  |_|      <------ (i,j,k = K): bottom layer
    |            | |
    |            |_|        <------ (i,j,k = 2): depth layer
    |            | 
    |____________|___ J     <------ (i,j,k = 1): reservoir's surface
    |    
    |
     I

%}

%----------------- POROSITY
disp('Rearranging data arrays...');

% TODO Change the counter 10 so as to generalize the method

B = [];    
for m = 1:J*K
    A = phi( (m-1)*10+1:10*m , : );    % loading matrix 10x6 from file
    A = A';                            % transpose to order
    A2 = reshape( A, [ 1 numel(A) ] ); % vector 1x60
    A2 = A2';                          % transpose to 60x1
    B = [ B; A2 ];                     % big column vector                           

end            
        
%{            
       [  [ j=1 ] [ j=1 ] . . . [ j=1 ]  ] ___ line 1*I 
       [  [ j=2 ] [ j=2 ]       [ j=2 ]  ] ___ line 2*I
       [     .       .    .        .     ]
       [     .       .      .      .     ]
 B11 = [     .       .        .    .     ]
       [  [ j=J ] [ j=J ] . . . [ j=J ]  ] ___ line J*I        
                                        
             |       |    . . .    |     
            k=1     k=2           k=K     
%}
B11 = reshape( B, [ I*J K ] );         

B22 = [];
for k = 1:K;
    for j = 1:J;
        B21 = B11( (j-1)*I+1:j*I, k );
        B22 = [ B22, B21 ];            
    end      
end

% filling porosity slices
PHI = zeros(I,J,K);
for m = 1:K
     B23 = B22( :, (m-1)*J+1:m*J );         
     PHI(:,:,m) = B23;         
end

disp('Porosity 3D array ready.');
%----------------- PERMEABILITY

%{
    In the case of porosity, the same scheme as B11 is 
    applied, except that we have 3 huge blocks now. 
    It looks like:
    
          [ B11x ]___ line 1*J*I
    B11 = [ B11y ]___ line 2*J*I
          [ B11z ]___ line 3*J*I
%}


% separating blocks
n = size(per,1)/3;
kx = per(     1:n   , : );
ky = per(   n+1:2*n , : );
kz = per( 2*n+1:end , : );

Bx = [];    
By = [];    
Bz = [];    
for m = 1:J*K
    Ax = kx( (m-1)*10+1:10*m , : );    % loading matrix 10x6 from file
    Ay = ky( (m-1)*10+1:10*m , : );    
    Az = kz( (m-1)*10+1:10*m , : );    
    Ax = Ax';                            % transpose to order
    Ay = Ay';                          
    Az = Az';                          
    A2x = reshape( Ax, [ 1 numel(Ax) ] ); % vector 1x60
    A2y = reshape( Ay, [ 1 numel(Ay) ] ); 
    A2z = reshape( Az, [ 1 numel(Az) ] ); 
    A2x = A2x';                          % transpose to 60x1
    A2y = A2y';                          
    A2z = A2z';                          
    Bx = [ Bx; A2x ];                     % big column vector                           
    By = [ By; A2y ];                    
    Bz = [ Bz; A2z ];                    

end            
        
B11x = reshape( Bx, [ I*J K ] );         
B11y = reshape( By, [ I*J K ] );         
B11z = reshape( Bz, [ I*J K ] );         

B22x = [];
B22y = [];
B22z = [];
for k = 1:K;
    for j = 1:J;
        B21x = B11x( (j-1)*I+1:j*I, k );
        B21y = B11y( (j-1)*I+1:j*I, k );
        B21z = B11z( (j-1)*I+1:j*I, k );
        B22x = [ B22x, B21x ];            
        B22y = [ B22y, B21y ];            
        B22z = [ B22z, B21z ];            
    end      
end

% filling permeability slices
KX = zeros(I,J,K);
KY = zeros(I,J,K);
KZ = zeros(I,J,K);
for m = 1:K
     B23x = B22x( :, (m-1)*J+1:m*J );         
     KX(:,:,m) = B23x;
     
     B23y = B22y( :, (m-1)*J+1:m*J );         
     KY(:,:,m) = B23y;
     
     B23z = B22z( :, (m-1)*J+1:m*J );         
     KZ(:,:,m) = B23z;
end

disp('Permeability 3D arrays ready.');

end

