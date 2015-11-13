%% mainVOIWellGraphData

% load DRT matrix
aux = load('../mat/DRT_Field.mat');
DRT = aux.DRT;
id = find(DRT(:) == -Inf); % eliminating -Inf
DRT(id) = 0.0;

% directory
matFiles = dir('../mat/DRT_VOI*.mat'); 
numfiles = length(matFiles);

% statistics
drtval = zeros(numfiles,1);
qntComps = zeros(numfiles,1); 
qntGreaterComps = zeros(numfiles,1);

% voxelCoords
bigNetworkCoords = cell(numfiles,1);

for k = 1:numfiles
    st = load( strcat('../mat/',matFiles(k).name) ); 
    val = st.VOISt.value;     
    drtval(k) = val;
    
    ncomp = length(st.VOISt.compVoxelCoords);
    qntComps(k) = ncomp;
    qntGreaterComps(k) = st.VOISt.compNNodes{1}; 
    
    bigNetworkCoords{k} = st.VOISt.compVoxelCoords{1};
    
    %for idcomp = 1:ncomp
    %    cvi = st.VOISt.compVoxelInds{idcomp};         
    %    plotVoxelGraphComp( DRT,cvi,val,idcomp,0.8);                
    %end   
    
    % biggest component 
    cvi = st.VOISt.compVoxelInds{1};         
    plotVoxelGraphComp( DRT,cvi,val,1,0.1,gray,'p', [-37.5,30],'pdf');                    
    
    
end

% statistics
M = [ drtval, qntComps, qntGreaterComps ];


