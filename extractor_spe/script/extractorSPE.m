%% extractorSPE.m - Data extractor for oil wells
%{     
    authors: Dr. Gustavo Peixoto de Oliveira
             Msc. Tatiana Araujo Simoes 
             Dr. Moises Dantas dos Santos
             Dr. Waldir Leite Roque
             @Federal University of Paraiba
    mail: gustavo.oliveira@ci.ufpb.br    
    date: Sep 21st, 2015        
            
    input: 
        .dat files for porosity and permeability according to SPE data 
    
    output:
        several data for analysis and visualization
    
%}
%% DEFAULTS

clear all; close all; clc; format long;
diary('extractorSPE.log');
diary on

setOptions;     % defaults
splshScreen;    % screen

%% LOAD FILES
phiname = '../dat/spe_phi.dat';
pername = '../dat/spe_perm.dat';
[phi,per] = loadFiles(phiname,pername);

%% GRID BOUNDS (SPE Project 2)
[I,J,K] = setGridBounds(60,220,85); % default


%% RUN OPTIONS

opt = input('----> Choose method for well extraction: [0] N random wells; [1] specific well: \n');

if ( opt == 0 || isempty(opt) ) % random
    
    N = input('----> Choose number of random wells to extract: \n');
    
    if N > I*J
        error( strcat(' > Maximum number of wells to extract is ', num2str(I*J),'!') ); 
    else        
        ia = randperm(I);   ia = ia(1:N);
        ja = randperm(J);   ja = ja(1:N);
        printWellTable( ia, ja, N );
    end
    
elseif opt == 1 % specific
    
    N = 1;
    
    ia = input('-------> Enter with I coordinate of the well: \n');
    ja = input('-------> Enter with J coordinate of the well: \n');    
              
    if ia > I || ja > J 
        error( strcat( 'Maximum I= ',num2str(I), ' Maximum J= ',num2str(J) ) );
    else
        printWellTable( ia, ja, N );
    end
    
else
    error('Option not valid.');
end

pltd = false;
pltda = input('----> Plot 3D reservoir? [0] no; [1] yes. \n');
if pltda == 1
    pltd = true;
end

pltr = false;
pltra = input('----> Plot and print well data to file? [0] no; [1] yes. \n');
if pltra == 1
    pltr = true;
end

vtk = false;
vtkd = input('----> Export data to VTK for 3D visualization? [0] no; [1] yes. \n');
if vtkd == 1;
    vtk = true;
end

csv = false;
csvd = input('----> Export well data to CSV? [0] no; [1] yes. \n');
if csvd == 1;
    csv = true;
end

wsati = 0.0; % irreducible water saturation distribution
%%%%%%%%%%%%       THIS WILL BE USED LATER!       %%%%%%%%%%%%%%%%%%%%%%%%%
%
% ws = input('----> Include irreducible water saturation distribution? [0] no; [1] yes. \n');
% if ws == 1;            
%     
%     wc = input('-------> Enter with default fraction values? [0] no; [1] yes. \n'); 
%     
%     if wc == 0
%         wsatimin1 = input('-------> Enter with minimum fraction for wsati for gas layer: \n'); 
%         wsatimax1 = input('-------> Enter with maximum fraction for wsati for gas layer: \n'); 
%         wsatimin2 = input('-------> Enter with minimum fraction for wsati for oil layer: \n'); 
%         wsatimax2 = input('-------> Enter with maximum fraction for wsati for oil layer: \n'); 
%         wsatimin3 = input('-------> Enter with minimum fraction for wsati for water layer: \n'); 
%         wsatimax3 = input('-------> Enter with maximum fraction for wsati for water layer: \n'); 
%            
%     else % paper by Chandra (2008)
%         wsatimin1 = 0.4;
%         wsatimax1 = 0.8;
%         wsatimin2 = 0.4;
%         wsatimax2 = 1.0;
%         wsatimin3 = 0.8;
%         wsatimax3 = 1.0;
%     end
%         
%     if wsatimin1 > wsatimax1 || wsatimin2 > wsatimax2 || wsatimin3 > wsatimax3
%         error('Minimum fraction > maximum fraction??');    
%     else
%         printWSatTable(wsatimin1,wsatimax1,wsatimin2,wsatimax2,wsatimin3,wsatimax3);
%     end               
%     
%     %{
%         Generating a normal random distribution for the 
%         irreducible water saturation profile along the 3 layers.
%     %}
%     
%     n1 = round(K/3);     
%     n3 = n1;
%     n2 = K - (n1+n3); % mid layer with preference.
%     dc = 0.2; % deviation constant (value of 30% generates wsati > 1 => imaginary DRT)
%     
%     % gas layer (top)
%     mean1 = 0.5*(wsatimin1 + wsatimax1);
%     stda1 = dc*(wsatimax1 - wsatimin1); 
%     wsati1 = mean1 + stda1*randn(n1,1);
% 
%     % oil layer (mid)
%     mean2 = 0.5*(wsatimin2 + wsatimax2);
%     stda2 = dc*( wsatimax2 - wsatimin2 );
%     wsati2 = mean2 + stda2*randn(n2,1);
%         
%     % water layer (bottom)
%     mean3 = 0.5*(wsatimin3 + wsatimax3);
%     stda3 = dc*( wsatimax3 - wsatimin3 );
%     wsati3 = mean3 + stda3*randn(n3,1);
%     
%     wsati = [ wsati1; wsati2; wsati3 ]; % profile
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bland = false;
blandd = input('----> Plot and export Bland-Altman diagrams to file? [0] no; [1] yes. \n');
if blandd == 1;
    bland = true;
end

hist = false;
histd = input('----> Plot and export histograms and regression analyses to file? [0] no; [1] yes. \n');
if histd == 1;
    hist = true;
end

svmat = false;
svmatd = input('----> Save .mat 3D arrays (porosity, permeability) to file? [0] no; [1] yes. \n');
if svmatd == 1;
    svmat = true;
end

pltdrt = false;
pltdrtd = input('----> Plot and export HFU locations to file? [0] no; [1] yes. \n');
if pltdrtd == 1;
    pltdrt = true;
end

disp('Options saved. Running extractor...');

%% MATRIX OPERATIONS

%{
    This strategy is ONLY applicable to the original .dat files 
    provided by the SPE Project 2 and follows the data arrangement.
    It will might be useful for other files as well, but further
    investigation is required. 

    What is done here: the original files are 'unrolled' 
    and rearranged into 3D arrays (I,J,K), thus allowing a 
    better accessibility of the information. In the end, the reservoir
    is a 3D structure with { PHI, KX, KY, KZ }(I,J,K) evaluated 
    in the fourth dimension. 

    RESERVOIR
    =========

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
    |____________|___ J     <------ (i,j,k = 1): reservoir surface
    |    
    |
     I

%}

%----------------- POROSITY
disp('Rearranging data arrays...');

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
             
%% WELL STRUCTURE

%{
    Since we may have N random wells (w) to extract at once, 
    the strategy is to store the desired information 
    as cell arrays. Later on, a big cell arrays 
    with its rows representing the wells
    and columns the computed properties is arranged: 


           { [phi] [kx] [ky] [kz] ... } <-- well 1  
           { [phi] [kx] [ky] [kz] ... } <-- well 2  
    wMat = {           .              }
           {           .              }
           {           .              } 
           { [phi] [kx] [ky] [kz] ... } <-- well N  
    
%}


disp('Creating cell arrays for wells...');
% storing cells 
wphi  = cell(1,N);
wkx   = cell(1,N);
wky   = cell(1,N);
wkz   = cell(1,N);
wphiz = cell(1,N);
wkn   = cell(1,N);
wRQI  = cell(1,N);
wFZI  = cell(1,N);
wDRT  = cell(1,N);
wMFZI = cell(1,N);
wMDRT = cell(1,N);
wFZIStar = cell(1,N);
wLogPhiz = cell(1,N);
wLogRQI = cell(1,N);
wMat  = cell(N,17);

disp('Sweeping wells...');
% sweeping wells
for i = 1:N    

    aux = reshape( PHI(ia(i),ja(i),:), [ K 1 ] );
    wphi{i} = aux; 
            
    wphiz{i} = wphi{i}./(1.0 - wphi{i}); % normalized porosity
    
    aux = reshape( KX(ia(i),ja(i),:), [ K 1 ] );
    wkx{i} = aux;    
    
    aux = reshape( KY(ia(i),ja(i),:), [ K 1 ] );
    wky{i} = aux;    
    
    aux = reshape( KZ(ia(i),ja(i),:), [ K 1 ] );
    wkz{i} = aux;    
    
    wkn{i} = sqrt( wkx{i}.^2 + wky{i}.^2 + wkz{i}.^2 ); % permeability norm
    
    wRQI{i} = 0.0314*sqrt( wkn{i}./wphi{i} );
    
    % AMAEFULE method (no water saturation)
    wFZI{i} = wRQI{i}./wphiz{i};            
    wDRT{i} = round( 2*log( wFZI{i} ) + 10.6 );
    
    % IG method (water saturation)
    wMFZI{i} = wRQI{i}./( wphiz{i}.*(1.0 - wsati).^(3.0/2.0) );    
    wMDRT{i} = round( 2*log( wMFZI{i} ) + 10.6 );
    
    % PJP method ( Energy Tech 2015 (3), 726-733 )
    wFZIStar{i} = wRQI{i};    

    % Logarithms (added to ease fit data plot)
    wLogPhiz{i} = log( wphiz{i} );
    wLogRQI{i} = log( wRQI{i} );
    
    % indices ( i = ia, j = ja, k )    
    wia = ia(i)*ones( length( wphi{i} ),1 );
    wja = ja(i)*ones( length( wphi{i} ),1 );
    wka = (1:K)';
    
    %{
        These operations eliminate all the entries 
        of the vectors above where porosity = 0 was found 
        to avoid indetermination.         
        Strategy:
            i)  find nonzero entries of porosity
            ii) shorten all the vectors by exclusion                    
    %}    
    
    % find
    wid = find( wphi{i} ~= 0 );     
    wid0 = find( wphi{i} == 0 );   
    saveExpurgatedIds( wid0, ia(i), ja(i) ); % export expurgated points
    
    % shorten
    wia = wia(wid);
    wja = wja(wid);
    wka = wka(wid);    
    wphi{i} = wphi{i}(wid); 
    wkx{i} = wkx{i}(wid); 
    wky{i} = wky{i}(wid); 
    wkz{i} = wkz{i}(wid); 
    wkn{i} = wkn{i}(wid); 
    wphiz{i} = wphiz{i}(wid); 
    wRQI{i} = wRQI{i}(wid);
    wFZI{i} = wFZI{i}(wid);
    wDRT{i} = wDRT{i}(wid);
    wMFZI{i} = wMFZI{i}(wid);
    wMDRT{i} = wMDRT{i}(wid);
    wFZIStar{i} = wFZIStar{i}(wid);
    wLogPhiz{i} = wLogPhiz{i}(wid);
    wLogRQI{i} = wLogRQI{i}(wid);
    
    % organize into cell as matrix whose rows contain data for each well
    wMat{i,1} = wia;
    wMat{i,2} = wja;
    wMat{i,3} = wka;                              
    wMat{i,4} = wphi{i};
    wMat{i,5} = wkx{i};
    wMat{i,6} = wky{i};
    wMat{i,7} = wkz{i};
    wMat{i,8} = wkn{i};
    wMat{i,9} = wphiz{i};
    wMat{i,10} = wRQI{i};
    wMat{i,11} = wFZI{i};
    wMat{i,12} = wDRT{i};
    wMat{i,13} = wMFZI{i};
    wMat{i,14} = wMDRT{i};
    wMat{i,15} = wFZIStar{i};
    wMat{i,16} = wLogPhiz{i};
    wMat{i,17} = wLogRQI{i};
    
    fprintf('Well %d: (%d,%d) done.\n',i,ia(i),ja(i) );
end


%% RESERVOIR PLOT

% mounting meshgrid for plotting     
i = 1:I;
j = 1:J;
k = K:-1:1; % contrary to show plot upward
[JJ,II,KK] = meshgrid(j,i,k);

if pltd    
    disp('Plotting 3D reservoirs...')
    figure 
    h = gca;
    set(h,'YDir','reverse'); % z-depth growing downward
    slice(JJ,II,KK,PHI,j,i,k)
        
    figure 
    h = gca;
    set(h,'YDir','reverse');
    slice(JJ,II,KK,KX,j,i,k)            
    
    figure 
    h = gca;
    set(h,'YDir','reverse');
    slice(JJ,II,KK,KY,j,i,k)            
    
    figure 
    h = gca;
    set(h,'YDir','reverse');
    slice(JJ,II,KK,KZ,j,i,k)
end

%----------------- DISPERSION GRAPHS 
%{

  depth
  z @surface (voxel 1)
    | o
    |     o
    |    
    |o   
    |  o 
    |       o
     ---------> phi, k(.)
         

%}

if pltr      
    
    %{ 
        Remark: fliplr(k) and 'Ydir reverse' are necessary to view 
        the dispersion data of the reservoir's with z = 0 at the
        free surface and z = zmax at the bottom.        
      
    %}
    dpx = 0.3; % graph padding
    dpy = 2;    
    figname = 'dispersion_I';
    
    i = 1;
    while i <= N     
        fprintf('Plotting dispersion graphs for well %d...\n',i)                  
        
        % --------------- porosity 
        figure                                
        subplot(2,2,1)              
        scatter(wphi{i},fliplr(wMat{i,3}),'fill','d','MarkerFaceColor','r')
        h = gca;
        set(h,'YGrid','on','YDir','reverse','YTick');
        xlim( [ (1.0 - dpx)*min(wphi{i}) (1.0 + dpx)*max(wphi{i}) ] )
        ylim( [ -dpy K+dpy ] );        
        title( strcat('Well (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
        xlabel('$ \phi $','interpreter','latex');
        ylabel('z'); 
        
        % --------------- porosity - x 
        subplot(2,2,2)                 
        scatter(wkx{i},fliplr(wMat{i,3}),'fill','d','MarkerFaceColor','g')
        h = gca;
        set(h,'YGrid','on','YDir','reverse','YTick');
        xlim( [ (1.0 - dpx)*min(wkx{i}) (1.0 + dpx)*max(wkx{i}) ] )
        ylim( [ -dpy K+dpy ] );
        title( strcat('Well (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
        xlabel('$ \kappa_x $','interpreter','latex');
        ylabel('z'); 
        
        % --------------- porosity - y
        subplot(2,2,3)                 
        scatter(wky{i},fliplr(wMat{i,3}),'fill','d','MarkerFaceColor','b')
        h = gca;
        set(h,'YGrid','on','YDir','reverse','YTick');
        xlim( [ (1.0 - dpx)*min(wky{i}) (1.0 + dpx)*max(wky{i}) ] )
        ylim( [ -dpy K+dpy ] );
        title( strcat('Well (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
        xlabel('$ \kappa_y $','interpreter','latex');
        ylabel('z'); 
        
        % --------------- porosity - z
        subplot(2,2,4)                        
        scatter(wkz{i},fliplr(wMat{i,3}),'fill','d','MarkerFaceColor','m')        
        h = gca;
        set(h,'YGrid','on','YDir','reverse','YTick');
        xlim( [ (1.0 - dpx)*min(wkz{i}) (1.0 + dpx)*max(wkz{i}) ] )
        ylim( [ -dpy K+dpy ] );        
        title( strcat('Well (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
        xlabel('$ \kappa_z $','interpreter','latex');
        ylabel('z'); 
                 
        % print to file
        print('-dpdf','-r0',fullfile( '../figs/', strcat(figname,num2str( ia(i) ),'_J', num2str( ja(i) ) ) ) );
        
        i = i+1; 
    end
end

%% VTK EXPORT

if vtk
    disp('Exporting to VTK...');
    savevtk_structured_spe(I,J,K,PHI,KX,KY,KZ,'../vtk/spe-phi-k-reservoir');                    
    
    i = 1;    
    while i <= N  
        I = ia(i); J = ja(i) ;
        savevtk_structured_well_spe( i,I,J,wphi{i},wkx{i},wky{i},wkz{i},...
                            strcat( '../vtk/spe-phi-k-well',num2str(i),'-',num2str(I),'_',num2str(J)) );                            
        i = i + 1;        
    end
end

%% CSV EXPORT
%{
 
Format of CSV file
==================

[ 'i' 'j' 'k' 'phi_e' 'kx' 'ky' 'kz' 'kn' 'phiz' 'RQI' 'FZI' 'DRT' 'MFZI' 'MDRT' 'FZI*' 'LogPhiz' 'LogRQI' ]; 

    i,j,k: triplet indices
    phi_e: effective porosity
       kx: permeability - x 
       ky: permeability - y
       kz: permeability - z
       kn: norm-2 of permeability 
     phiz: normalized porosity
      RQI: Reservoir Quality Index
      FZI: Flow Zone Indicator
      DRT: Discrete Rock Type
   M(...): Modified (...)
     FZI*: method PJP
  LogPhiz: log(phiz)
   LogRQI: log(RQI)
    
%}
if csv        
    csvname = 'well_I';        
    i = 1;
    while i <= N    
        aux = [ wMat{i,1} wMat{i,2} wMat{i,3} wMat{i,4} ... 
                wMat{i,5} wMat{i,6} wMat{i,7} wMat{i,8} ... 
                wMat{i,9} wMat{i,10} wMat{i,11} wMat{i,12} ...
                wMat{i,13} wMat{i,14} wMat{i,15} wMat{i,16} wMat{i,17} ];
        fprintf('Exporting CSV file for well %d... \n',i);  
        fname = fullfile( '../csv/', strcat(csvname,num2str( ia(i) ),'_J', num2str( ja(i) ),'.csv' ) );
        csvwrite(fname, aux);
        
        % TODO
        %head = cellstr( ['i     '; 'j     '; 'k     '; 'phi_e '; 'kx    '; 'ky    '; 'kz    '; 'kn    '; 'phiz  '; 'RQI   '; 'FZI   '; 'DRT   '; 'MFZI  '; 'MDRT  '; 'FZI*  '; 'logPhiz'; 'logRQI'] ); 
        %head = head';
        %csvwrite_with_headers(fname,head,aux); % complaining about cellstr.
        i = i + 1;
    end
end


%% Saving PHI,KX=KY,KZ to file for posterior use

if svmat
    disp('Saving .mat files...');
    save('../mat/PHI.mat','PHI');
    save('../mat/KX.mat','KX');
    save('../mat/KY.mat','KY');
    save('../mat/KZ.mat','KZ');
end


%% Bland-Altman 

if bland && ~isempty(wsati) % if wsat is empty, there's nothing to compare
    
    dif = cell(N,1);
    med = cell(N,1);
    medDif = zeros(N,1);
    stdDif = zeros(N,1);
    sup = zeros(N,1);
    inf = zeros(N,1);

    baname = 'Bland-Altman_I';

    for i = 1:N
        fprintf('Saving Bland-Altman diagram for well %d... \n',i);  
        
        dif{i} = wMat{i,11} - wMat{i,13};    
        med{i} = 0.5*( wMat{i,11} + wMat{i,13} );
        medDif(i) = mean( dif{i} );
        stdDif(i) =  std( dif{i} );

        sup(i) = medDif(i) + 1.96*stdDif(i);
        inf(i) = medDif(i) - 1.96*stdDif(i);

        figure 
        hold on    
        scatter(med{i},dif{i});
        xlim( [ -0.5 0.1*max(med{i})+0.5] );
        ylim( [ inf(i)-0.5 sup(i)+0.5] );
        line([0 0.1*max(med{i}) ], [medDif(i) medDif(i)],'Marker','.','LineStyle','-')
        line([0 0.1*max(med{i}) ], [sup(i) sup(i)],'Marker','.','LineStyle','--')
        line([0 0.1*max(med{i}) ], [inf(i) inf(i)],'Marker','.','LineStyle','--')
        title( strcat('Bland-Altman (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );        
        xlabel('ave');
        ylabel(' $ FZI_{AM-IG} $','interpreter','latex');  

        % print to file
        print('-dpdf','-r0',fullfile( '../figs/', strcat(baname,num2str( ia(i) ),'_J', num2str( ja(i) ) ) ) );

    end
    
end

%% Histogram, Regression Analysis

if hist

    hname = 'HistogramDRTs';
    regname = 'Regression_I';
    nfreq = 5; % number of DRT frequencies to analyze in regression 
    seps = 0.05; % slope tolerance for regression fit line (1.0-seps,1.0+seps)

    for i = 1:N                        
        tbh = tabulate( wMat{i,12} ); % returns 'histogram' data:
        idt = find( tbh(:,2) ~= 0 ); % excludes 0 frequencies

        drt = tbh(:,1); 
        drt = drt(idt); % DRT value ~= 0
        frq = tbh(:,2);
        frq = frq(idt); % statistical frequencies

        figure    
        stem(drt,frq,'fill','k');
        grid('on');    
        xlim([ min(drt)-0.5 max(drt)+0.5 ]);
        title( strcat('Histogram - DRTs (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );            
        xlabel(' $ DRT $','interpreter','latex');  
        ylabel('frequency');
        print('-dpdf','-r0',fullfile( '../figs/', ...
                       strcat(hname,num2str( ia(i) ),'_J', ...
                                      num2str( ja(i) ) ) ) ); 

        % create cell for best DRTs found                         
        if i == 1
            dataDRT = cell( length(drt), 6 );        
        end

        drtgood = [];  % best DRTs
        for j = 1:length(drt)
            idk = find( wMat{i,12} == drt(j) );
            dataDRT{j,1} = log( wMat{i,9}(idk) );  % phiz
            dataDRT{j,2} = log( wMat{i,10}(idk) ); % RQI
            dataDRT{j,3} = wMat{i,3}(idk);         % z coordinates        

            [rj, mj, bj ] = regression( dataDRT{j,1}, dataDRT{j,2},'one' );

            dataDRT{j,4} = rj; % regression value (Pearson correlation)                    
            dataDRT{j,5} = mj; % slope              
            dataDRT{j,6} = bj; % offset        

            % plot only DRT with a good frequency and slopes close to 1.0
            % by slpeps
            if size( dataDRT{j,1}, 1 ) >= nfreq && mj >= 1.0 - seps && mj <= 1.0 + seps          
                
                drtgood = [ drtgood; drt(j) ];                
                
                figure             
                PLR = plotregression( dataDRT{j,1}, dataDRT{j,2}, ...
                    strcat('(',num2str( ia(i) ),    ...
                                   ',',num2str( ja(i) ),')',...
                                ' DRT=',num2str( drt(j) ),' tol_s=',num2str(seps) ) );
                
                setLRPlot(PLR); % graph appearance
                
                % csv file [ logphiz, logRQI, depth ] per DRT                
                fprintf('Exporting CSV file of logs/depth data for DRT %d... \n', drt(j) );  
                fname = fullfile( '../csv/', strcat(regname,num2str( ia(i) ),'_J', num2str( ja(i) ),'_DRT_',num2str( drt(j) ),'_LogsDepth','.csv' ) );
                auxmat = [ dataDRT{j,1} dataDRT{j,2} dataDRT{j,3} ];
                csvwrite(fname, auxmat );                                

                % csv file: fit data
                fprintf('Exporting CSV file of regression fit data for DRT %d... \n', drt(j) );  
                fname = fullfile( '../csv/', strcat(regname,num2str( ia(i) ),'_J', num2str( ja(i) ),'_DRT_',num2str( drt(j) ),'_FitData','.csv' ) );
                aux = [ dataDRT{j,4} dataDRT{j,5} dataDRT{j,6} ];
                csvwrite(fname, aux );
                
                % print to file    
                print('-dpdf','-r0',fullfile( '../figs/', ...
                       strcat(regname,num2str( ia(i) ),'_J', ...
                                      num2str( ja(i) ), ...
                                      '_DRT_',num2str( drt(j) ) ) ) );                
                
            end

        end   
        
        % TODO (OK for a unique well; not optimized for all)
        % execCSVCondenser; % merging .csv files of the best DRTs into one
        
        
        if pltdrt
            % Selecting the best DRTs for plotting
            depths = cell(1,length(drtgood));
            for k = 1:length(drtgood)
                id = drtgood(k);
                idd = find( drt(:) == id );
                depths{1,k} = dataDRT{ idd,3 };
            end

            % plot the best HFU locations per well
            plotHFULoc( drtgood,depths,ia(i),ja(i) ); 
        end
                

    end

end

%% Final output
disp('Closing figures...')
close all;
disp('----> E X E C U T I O N    H A S    F I N I S H E D.');

diary off