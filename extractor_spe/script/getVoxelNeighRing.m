function [FN,VN,FZ,VZ] = getVoxelNeighRing(ic,jc,kc,P,F)
%{ 
    getVoxelNeighRing - gets the voxel neighborhood 
    VN(ic,jc,kc;P) by a ring of P of the 3D array
    F. The neighborhood is a cube growing with 
   (2*P+1)^3 voxels
    
    input: 
      ic,jc,kc: central voxel coordinate
             P: voxel ring radius
             F: scalar field values (3D array)
    output: 
        VN: voxel global indices
        FN: scalar field evaluated in the neighborhood
        FZ: voxel global indices whose entries are zero (array)
        VZ: voxel local indices whose entries are zero (neighborhood)

    SCHEME
    ======

    - Layer k = kc
        
                 jc+P
                 .
                 .
                 .
                          ____
                         |    |
                         |    |
                ____ ____|____|
               |    |    |    |
   ic-P <---   | c  |c+1 |c+2 | ---> ic+P
               |____|____|____|
               |    |    |    |
               | c-1|    |    |
           ____|____|____|____| 
          |    |    |    |    |
          |    | c-2|    |    |  
          |____|____|____|____|
                
                 .               for k = kc-P:kc+P
                 .
                 .
                 jc - P                                       
     
  Gustavo Peixoto, 8 Oct 2015, @UFPB

%}

if P <= 0, error('P-ring must be > 0'); end

[I,J,K] = size(F);

vx = (2*P + 1); % voxel radius

%{ 
    Check for bounds:
    (ic,jc,kc) must be chosen in such a way 
    that the P-ring doesn't leak out the 
    boundaries    
%}
if ( ic < vx     || jc < vx     || kc < vx     || ... 
     ic > (I-vx) || jc > (J-vx) || kc > (K-vx)  )
    error('Voxel neighborhood ring exceeds the boundaries. Change central coordinate.')
end

nvx = vx^3;  % total number of voxels
vx3 = nthroot(nvx,3);

% voxel neighborhood (i,j,k) triplet matrix 
VN = [];

% sweeping voxel neighborhood per layer
for KP = kc-P:kc+P
    for IP = ic-P:ic+P
        for JP = jc-P:jc+P        
            VN = [ VN; IP JP KP ];
        end
    end
end

% evaluating scalar field
FN = zeros(vx3,vx3,vx3);
FZ = []; VZ = []; 
count = 1;
for k = 1:vx3
    for i = 1:vx3
        for j = 1:vx3    
            
            % scalar field
            FN(i,j,k) = F( VN(count,1), VN(count,2), VN(count,3) );
            
            % store voxel indices whose scalar value is 0.0
            if FN(i,j,k) == 0.0 
                % global
                FZ = [FZ; [ VN(count,1) VN(count,2) VN(count,3) ] ];
                % local
                VZ = [VZ; [ i j k ] ];
            end            
            count = count + 1;
        end
    end
end

end

           