function [ FN, VC, VN ] = getVoxelNeigh( ic, jc, kc, F )

[I, J, K] = size(F);

t = I - 3;
u = J - 3;
v = K - 3;

if ( ic < 3 || jc < 3 || kc < 3 || ... 
     ic > t || jc > u || kc > v  )
    error('Voxel neighborhood not applicable near the boundaries.')
end

%%%%%%%%% scalar field array

FN = zeros(3,3,3); 

%-------------layer k = kc-1

FN(1,1,1) = F( ic-1, jc-1, kc-1 );
FN(1,2,1) = F( ic-1, jc  , kc-1 );
FN(1,3,1) = F( ic-1, jc+1, kc-1 );

FN(2,1,1) = F( ic  , jc-1, kc-1 );
FN(2,2,1) = F( ic  , jc  , kc-1 );
FN(2,3,1) = F( ic  , jc+1, kc-1 );

FN(3,1,1) = F( ic+1, jc-1, kc-1 );
FN(3,2,1) = F( ic+1, jc  , kc-1 );
FN(3,3,1) = F( ic+1, jc+1, kc-1 );


%-------------layer k = kc

FN(1,1,2) = F( ic-1, jc-1, kc );
FN(1,2,2) = F( ic-1, jc  , kc );
FN(1,3,2) = F( ic-1, jc+1, kc );

FN(2,1,2) = F( ic  , jc-1, kc );
FN(2,2,2) = F( ic  , jc  , kc );
FN(2,3,2) = F( ic  , jc+1, kc );

FN(3,1,2) = F( ic+1, jc-1, kc );
FN(3,2,2) = F( ic+1, jc  , kc );
FN(3,3,2) = F( ic+1, jc+1, kc );


%-------------layer k = kc+1

FN(1,1,3) = F( ic-1, jc-1, kc+1 );
FN(1,2,3) = F( ic-1, jc  , kc+1 );
FN(1,3,3) = F( ic-1, jc+1, kc+1 );

FN(2,1,3) = F( ic  , jc-1, kc+1 );
FN(2,2,3) = F( ic  , jc  , kc+1 );
FN(2,3,3) = F( ic  , jc+1, kc+1 );

FN(3,1,3) = F( ic+1, jc-1, kc+1 );
FN(3,2,3) = F( ic+1, jc  , kc+1 );
FN(3,3,3) = F( ic+1, jc+1, kc+1 );


%%%%%%%%% voxel indices array

VN = zeros(27,3); 

%{
    Layer k = kc-1
                        j
    [ 1 2 3 ]        --> 
    [ 4 5 6 ]       |   
    [ 7 8 9 ]        i
%}
VN(1,1:3) = [ic-1, jc-1, kc-1 ];
VN(2,1:3) = [ic-1, jc  , kc-1 ];
VN(3,1:3) = [ic-1, jc+1, kc-1 ];
VN(4,1:3) = [ic  , jc-1, kc-1 ];
VN(5,1:3) = [ic  , jc  , kc-1 ];
VN(6,1:3) = [ic  , jc+1, kc-1 ];
VN(7,1:3) = [ic+1, jc-1, kc-1 ];
VN(8,1:3) = [ic+1, jc  , kc-1 ];
VN(9,1:3) = [ic+1, jc+1, kc-1 ];

%{
    Layer k = kc

    [ 10 11 12 ]        
    [ 13 14 15 ] 
    [ 16 17 18 ]
%}

VN(10,1:3) = [ic-1, jc-1, kc ];
VN(11,1:3) = [ic-1, jc  , kc ];
VN(12,1:3) = [ic-1, jc+1, kc ];
VN(13,1:3) = [ic  , jc-1, kc ];
VN(14,1:3) = [ic  , jc  , kc ];
VN(15,1:3) = [ic  , jc+1, kc ];
VN(16,1:3) = [ic+1, jc-1, kc ];
VN(17,1:3) = [ic+1, jc  , kc ];
VN(18,1:3) = [ic+1, jc+1, kc ];

%{
    Layer k = kc+1
                           j
    [ 19 20 21 ]        --> 
    [ 22 23 24 ]       |   
    [ 25 26 27 ]        i
%}
VN(19,1:3) = [ic-1, jc-1, kc+1 ];
VN(20,1:3) = [ic-1, jc  , kc+1 ];
VN(21,1:3) = [ic-1, jc+1, kc+1 ];
VN(22,1:3) = [ic  , jc-1, kc+1 ];
VN(23,1:3) = [ic  , jc  , kc+1 ];
VN(24,1:3) = [ic  , jc+1, kc+1 ];
VN(25,1:3) = [ic+1, jc-1, kc+1 ];
VN(26,1:3) = [ic+1, jc  , kc+1 ];
VN(27,1:3) = [ic+1, jc+1, kc+1 ];

% central voxel index
VC = VN(14,:);

end