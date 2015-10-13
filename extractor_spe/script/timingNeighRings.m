function [ T, ttot ] = timingNeighRings( PHI,KX,KY,KZ,plt,fig,varargin )
%{
    timingNeighRings - Computes the time spent to perform the 
                       operations required to get the DRT values 
                       by considering all the allowed ring radii 
                       for the porosity, permeability 3D arrays

    input: 
       PHI,KX,KY,KZ: 3D array matrices for porosity, permeability
                plt: boolean to enable plot
                fig: boolean to enable figure print
    output: 
                  T: vector of times per ring
               ttot: total time spent
%}


if nargin < 4 
    error('Missing arguments'); 
elseif nargin == 4
    plt = false;
    fig = false;
elseif nargin == 5
    fig = false;
end

[I, J, K] = size(PHI);

[PMax,IC,JC,KC] = getMaxRing(I,J,K); % maximum ring

disp('Computing time...');

T = []; 
for P=1:PMax    
    tstart = tic; % marking time for these ops per ring

    [PHIV,~,~,~] = getVoxelNeighRing(IC,JC,KC,P,PHI);
    [KXV,~,~,~] = getVoxelNeighRing(IC,JC,KC,P,KX);
    [KYV,~,~,~] = getVoxelNeighRing(IC,JC,KC,P,KY);
    [KZV,~,~,~] = getVoxelNeighRing(IC,JC,KC,P,KZ);

    KXVN = sqrt( KXV.^2 + KYV.^2 + KZV.^2 );
    PHIVZ = PHIV./(1.0 - PHIV);
    RQIV = 0.0314*sqrt( KXVN./PHIV );
    FZIV = RQIV./PHIVZ;
    round( 2*log( FZIV ) + 10.6 ); % DRT

    telapsed = toc(tstart);
    T = [ T; telapsed ]; % time per iteration

end

ttot = cumsum(T); ttot = ttot(end); % total time
rings = 1:PMax; % ring radii vector

% dat file
fname = '../dat/subdomain/timesPerRing.dat';
disp( strcat('Time file saved to: ', fname) );
dlmwrite(fname,[ rings' T ], '\t');  

% plot 
if plt && ~fig
    figure
    bar(rings',T);
    title('Computational time');
    xlabel('rings');
    ylabel('time (s)');
end

% figure
if plt && fig
    figure
    bar(rings',T);
    title('Computational time');
    xlabel('rings');
    ylabel('time (s)');
    print('-depsc2','-r0','../figs/subdomain/timeRingGraph.eps');    
end


end

