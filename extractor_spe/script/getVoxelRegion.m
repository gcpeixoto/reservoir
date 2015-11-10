function [FN,VN,IND] = getVoxelRegion(ic,jc,kc,P,F)
%   GETVOXELREGION gets the voxel neighbourhood - 
%         region of interest - VN(ic,jc,kc;P) by a ring 
%         of P = [Pi, Pj, Pk] of the 3D array F. The neighborhood is a 
%         rectangle of dimensions Pi x Pj x Pk, with Pi,Pj,Pk
%         the ring for each direction.
%     
%     input: 
%       ic,jc,kc: central voxel coordinate
%              P: voxel ring radius array 
%              F: scalar field values (3D array)
%     output: 
%
%         FN: scalar field evaluated in the voxel
%         VN: voxel global coordinates
%        IND: voxel global linear indices
%
% 
%     SCHEME
%     ======
% 
%     - Layer k = kc
%         
%                  jc+Pj
%                  .
%                  .
%                  .
%                           ____
%                          |    |
%                          |    |
%                 ____ ____|____|
%                |    |    |    |
%   ic-Pi <---   | c  |c+1 |c+2 | ---> ic+Pi
%                |____|____|____|
%                |    |    |    |
%                | c-1|    |    |
%                |____|____|____|            
%                 
%                  .               for k = kc-Pk:kc+Pk
%                  .
%                  .
%                  jc-Pj                                       
%      
%   Dr. Gustavo Peixoto, 4 Nov 2015, @UFPB

% checking
if ~all(P) || any(P < 0)
    error('Radius must be > 0.');
end
    
% dealing with ring dimension
if length(P) == 1 && P > 0
    Pi = P;
    Pj = 0;
    Pk = 0;
    
elseif length(P) == 2
    Pi = P(1);
    Pj = P(2);
    Pk = 0;
    
elseif length(P) == 3
    Pi = P(1);
    Pj = P(2);
    Pk = P(3);
    
else
    error('Check P.');
end

% setting limits of the region of interest
Pia = Pi; Pib = Pi;    
Pja = Pj; Pjb = Pj;
Pka = Pk; Pkb = Pk;

if ic == 1
    Pia = 0; Pib = Pi;    
elseif ic == size(F,1)
    Pia = Pi; Pib = 0;
end

if jc == 1        
    Pja = 0; Pjb = Pj;    
elseif jc == size(F,2)    
    Pja = Pj; Pjb = 0;    
end

if kc == 1    
    Pka = 0; Pkb = Pk;
elseif kc == size(F,3)    
    Pka = Pk; Pkb = 0;
end

% voxel neighborhood (i,j,k) triplet matrix 
VN = [];

% sweeping voxel neighborhood per layer
IP = ic-Pia:ic+Pib;
JP = jc-Pja:jc+Pjb; 
KP = kc-Pka:kc+Pkb;

for kp = KP(1):KP(end)
    for ip = IP(1):IP(end)
        for jp = JP(1):JP(end)
            VN = [ VN; ip jp kp ];
        end
    end
end

% evaluating scalar field
FN = zeros( size(VN,1), 1 );
 
% scalar field
for i = 1:size(VN,1);
    FN(i) = F( VN(i,1), VN(i,2), VN(i,3) );                
end

% voxel global linear indices
IND = zeros( size(VN,1),1 );
for i = 1:size(VN,1)        
    IND(i) = sub2ind( size(F),VN(i,1),VN(i,2),VN(i,3) );
end

end

           