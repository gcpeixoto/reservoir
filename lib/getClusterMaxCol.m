function [ ilims,jlims,klims,leni,lenj,lenk,Ws,nvcol,cvcol] = getClusterMaxCol( cvc, DRT, val )
% GETCLUSTERMAXCOL gets data of voxels at highly filled columns in a 
%                  given cluster.
%                   
%   input:
%           cvc: cluster voxel coordinates
%           DRT: 3D DRT field array
%           val: DRT value of interest
%
%   output:
%           ilims,jlims,klims,
%           leni,lenj,lenk: cluster bounding box limits and lengths
%           Ws: surface coordinates of the wells where the max cols were
%               found
%           nvcol: number of voxels at the highly filled col
%           cvcol: voxel coords at the col

% bounding box limits
im=min(cvc(:,1)); iM=max(cvc(:,1)); ilims = [im,iM];
jm=min(cvc(:,2)); jM=max(cvc(:,2)); jlims = [jm,jM];
km=min(cvc(:,3)); kM=max(cvc(:,3)); klims = [km,kM];

% lenghts
leni = length(im:iM); lenj = length(jm:jM); lenk = length(km:kM);

%----- sweeping cluster box

a = 1; % well counter
LEND = [];
for i=im:iM
    for j=jm:jM           
        ID = [];                        
        for k=km:kM                     % depths
            base = [i,j,k];             
            II = strmatch(base,cvc);    % voxel belongs to cluster?
            if ~isempty(II)             % if yes, stores the depth
                ID = [ ID; II];                              
            end
        end
        qts{a} = ID;                    % number of voxels at cluster per well
        LEND = [ LEND; [i,j,length(qts{a})] ];
        a = a+1;
    end
end    

row = find( LEND(:,3) == max(LEND(:,3) ) );     % wells of maximum filling
Ws = [LEND(row,1),LEND(row,2)];                 % well surface coordinates

lws = size(Ws,1);                % number of wells              
nvcol = cell(lws,1);            % number of voxels per highly filled wells 
cvcol = cell(lws,1);             % voxel coords per highly filled wells 
kk = km:kM;

% sweeping highly filled wells inside cluster
for nw = 1:lws
    u = DRT( Ws(nw,1), Ws(nw,2), kk ) == val; 
    zz = find( u == 1 );    
    nvcol{nw} = length(zz);    
    is = Ws(nw,1)*ones(length(zz),1); 
    js = Ws(nw,2)*ones(length(zz),1); 
    ks = kk(zz)';
    cvcol{nw} = [is,js,ks];       
end


