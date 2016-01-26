function [ P ] = assemble3DPressure( press,I,J,K )
%     ASSEMBLE3DPRESSURE assembles 3D array for pressure
% 
%     input:
%         press: reservoir's porosity matrix in CMG output%         
%         I,J,K: reservoir's length, width, depth (voxel-based)
% 
%     output:
%         P: 3D arrays storing the pressure field data.
% 
%     - REMARK
%     
%     This strategy is, at first glance, ONLY applicable to 
%     .dat (.txt) files that follow the output format of CMG software. 
%     The current version is tested for the SPE Project 2 
%     data arrangement. It will might be useful for other 
%     files as well, but further investigation is required. 
% 
%     What is done here: the original files are 'unrolled' 
%     and rearranged into 3D arrays (I,J,K), thus allowing a 
%     better accessibility of the information. In the end, the reservoir
%     is a 3D structure with P(I,J,K) evaluated 
%     in the fourth dimension. 
% 
%     RESERVOIR SCHEME
%     ================
% 
%           K
%          /
%         /____________
%        /|            |
%       /_|__________  |
%      /|            | |
%     /_|__________  |_|      <------ (i,j,k = K): bottom layer
%     |            | |
%     |            |_|        <------ (i,j,k = 2): depth layer
%     |            | 
%     |____________|___ J     <------ (i,j,k = 1): reservoir's surface
%     |    
%     |
%      I
%
% -----------------------
% Dr. Gustavo Peixoto



%----------------- PRESSURE
disp('Rearranging data arrays...');

% TODO Change the counter 10 so as to generalize the method

B = [];    
for m = 1:J*K
    A = press( (m-1)*10+1:10*m , : );    % loading matrix 10x6 from file
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
P = zeros(I,J,K);
for m = 1:K
     B23 = B22( :, (m-1)*J+1:m*J );         
     P(:,:,m) = B23;         
end

disp('Pressure 3D array ready.');

end

