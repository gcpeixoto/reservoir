%% mainExtractorSPE.m - Data extractor for oil wells
%      
%     authors: Dr. Gustavo Peixoto de Oliveira
%              Dr. Waldir Leite Roque
%              @Federal University of Paraiba
%     mail: gustavo.oliveira@ci.ufpb.br    
%     date: Sep 21st, 2015        
%             
%     input: 
%         .dat files for porosity and permeability according to SPE model data 
%     
%     output:
%         several data for analysis and visualization
%     

%% DEFAULTS 

clear all; close all; clc;

% classes
dm = SPEDirManager;
dm.activateLog(mfilename);

d = SPEDisplay;
d.printSplScreen(mfilename); 
d.printings(d.author1,d.author2,d.inst,d.progStat{1});
d.setOptions;                
d.extractorSPEwarning;       

%% ADVICE ABOUT CALLING PRINT/PLOT FUNCTIONS 
%
% 'print' and 'plot' functions are time-consuming and were 
% switched off by default. In case of print/plot something to file
% or only for visualization purposes, switch on the required flags

% print flags 
pflag_dispersion = false;      % well's dispersion plots
pflag_baltman    = false;      % Bland-Altman plots  
pflag_regression = false;      % linear regression plots
pflag_histDRT    = false;      % DRT histogram plots

% plot flags 
pltflag_reservoir3D = false;   % 3D view of reservoir
pltflag_wellDisp    = false;   % well's dispersion plots (by depth)
pltflag_vtk         = false;   % export data to VTK
pltflag_HFULoc      = true;   % DRT map overview (HFU locations)
pltflag_regression  = false;   % regression fit-line plots
pltflag_histDRT     = false;   % histogram of DRT distribution


%% PATHS TO REQUIRED INPUT FILES (POROSITY AND PERMEABILITY)

% .mat files

phiname = '../mat/PHI.mat';
kxname  = '../mat/KX.mat';
kyname  = '../mat/KY.mat';
kzname  = '../mat/KZ.mat';

%. dat files
phidat = setFile('../dat/spe_phi.dat');
perdat = setFile('../dat/spe_perm.dat');    

%% GRID BOUNDS
[I,J,K] = setGridBounds(60,220,85);  % SPE 2 model default. DO NOT CHANGE!
    
%% RUN OPTIONS

% rerun / reload
inrr = input(d.extSPERerun);

if inrr == 1             % rerun       
                     
    d.dispMsg(d.extSPERebuild);        
    phi = load(phidat,'-ascii');
    per = load(perdat,'-ascii');
    [ PHI, KX, KY, KZ ] = assemble3DArrays( phi, per, I, J, K );
    
    svmat = true; % enables .mat saving
         
elseif inrr == 0         % reload    
             
    d.dispMsg(d.extSPEReload);                            
    load(phiname,'PHI');    
    load(kxname,'KX');
    load(kyname,'KY');
    load(kzname,'KZ');    
    
    svmat = false; % disables .mat saving
                      
else
    d.dispNotValid;
end

% method for well extraction: random / specific 
opt = input(d.extSPEWExt);

if ( opt == 0 || isempty(opt) ) % random    
    
    N = input(d.extSPEChooseN);    
    
    if N > I*J
        error( strcat(' > Maximum number of wells to extract is ', num2str(I*J),'!') ); 
    else        
        ia = randperm(I);   ia = ia(1:N);
        ja = randperm(J);   ja = ja(1:N);
        printWellTable( ia, ja, N );
    end
    
elseif opt == 1 % specific
    
    N = 1;
    
    ia = input(d.extSPEWellI);
    ja = input(d.extSPEWellJ);    
              
    if ia > I || ja > J 
        error( strcat( 'Maximum I= ',num2str(I), ' Maximum J= ',num2str(J) ) );
    else
        printWellTable( ia, ja, N );
    end
    
else
    d.dispNotValid;
end

% well data (gross csv)
csv = false;
csvd = input(d.extSPEGrossCSV);
if csvd == 1;
    csv = true;
end

% Regression analysis
hist = false;
histd = input(d.extSPELinRegr);
if histd == 1;
    hist = true;
end

%%%%%%%%%%%%%%%%% BLAND-ALTMAN PLOTS (UNDER DEVELOPMENT) %%%%%%%%%%%%%%%%%
bland = false;
%blandd = input('----> Compute Bland-Altman plots? [0] no; [1] yes. \n');
%if blandd == 1;
%    bland = true;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%    CONNATE WATER (UNDER DEVELOPMENT     %%%%%%%%%%%%%%%%%%%%%
wsati = 0.0; % irreducible water saturation distribution
%
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

disp('Options saved. Running extractor...');
             
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


%% 3D RESERVOIR PLOT

if pltflag_reservoir3D == true           
    plot3DReservoir(I,J,K,PHI); % visualize the 3D field (default: PHI)
end

%% DISPERSIONS BY DEPTH 

if pltflag_wellDisp == true      
            
    i = 1;
    while i <= N     
        fprintf('Plotting dispersions for well %d...\n',i)                  
                
        figure                                
        subplot(2,2,1)              
        plotDispersion( wphi{i},wMat{i,3},'phi',ia,ja,i,'d','r');                        
        subplot(2,2,2)                 
        plotDispersion( wkx{i},wMat{i,3},'kx',ia,ja,i,'d','g');                                      
        subplot(2,2,3)                 
        plotDispersion( wky{i},wMat{i,3},'ky',ia,ja,i,'d','b');                                        
        subplot(2,2,4)                        
        plotDispersion( wkz{i},wMat{i,3},'kz',ia,ja,i,'d','m');                        
                 
        % print to file
        if pflag_dispersion == true
            print('-dpdf','-r0',fullfile( '../figs/', strcat('Dispersion_I',num2str( ia(i) ),'_J', num2str( ja(i) ) ) ) );
        end
        
        i = i+1; 
    end
end

%% VTK EXPORT

if pltflag_vtk == true
    disp('Exporting to VTK...');
    savevtk_structured_spe(I,J,K,PHI,KX,KY,KZ,'../vtk/spe-phi-k-reservoir');                                
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
    
    % file header used in the loop
    head = {'i,'; 'j,'; 'k,'; 'phi_e,';            ...
           'kx,'; 'ky,'; 'kz,'; 'kn,'; 'phiz,';    ...
           'RQI,'; 'FZI,'; 'DRT,';'MFZI,';'MDRT,';'FZI*,'; ...
           'LogPhiz,'; 'LogRQI'}; 
    head = head';
    txt=sprintf('%s\t',head{:});
    txt(end)='';
    
    csvname = 'well_I';        
    i = 1;
    while i <= N    
        aux = [ wMat{i,1} wMat{i,2} wMat{i,3} wMat{i,4} ... 
                wMat{i,5} wMat{i,6} wMat{i,7} wMat{i,8} ... 
                wMat{i,9} wMat{i,10} wMat{i,11} wMat{i,12} ...
                wMat{i,13} wMat{i,14} wMat{i,15} wMat{i,16} wMat{i,17} ];
        
        fprintf('Exporting CSV file for well %d... \n',i);  
        fname = fullfile( '../csv/', strcat(csvname,num2str( ia(i) ),'_J', num2str( ja(i) ),'.csv' ) );        
        dlmwrite(fname,txt,'');
        dlmwrite(fname,aux,'-append');
                
        i = i + 1;
    end
end


%% Saving PHI,KX=KY,KZ to file for posterior use

if svmat == true
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
        fprintf('Saving Bland-Altman plot for well %d... \n',i);  
        
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
        if pflag_baltman == true
            print('-dpdf','-r0',fullfile( '../figs/', strcat(baname,num2str( ia(i) ),'_J', num2str( ja(i) ) ) ) );
        end
    end
    
end

%% Histogram, Regression Analysis and related CSV files

if hist
        
    hname = 'HistogramDRTs';
    regname = 'Regression_I';
    nfreq = 3; % number of DRT frequencies to analyze in regression 
    seps = 0.05; % slope tolerance for regression fit line (1.0-seps,1.0+seps)

    for i = 1:N                        
        tbh = tabulate( wMat{i,12} ); % returns 'histogram' data:
        idt = find( tbh(:,2) ~= 0 ); % excludes 0 frequencies

        drt = tbh(:,1); 
        drt = drt(idt); % DRT value ~= 0
        frq = tbh(:,2);
        frq = frq(idt); % statistical frequencies

        if pltflag_histDRT == true
            figure    
            bar(drt,frq,'k');        
            xlim([ min(drt)-0.5 max(drt)+0.5 ]);
            title( strcat('Histogram - DRTs (',num2str( ia(i) ),',',num2str( ja(i) ),')' ) );            
            xlabel(' $ DRT $','interpreter','latex');  
            ylabel('frequency');
        end
        
        if pflag_histDRT == true
            print('-dpdf','-r0',fullfile( '../figs/', ...
                       strcat(hname,num2str( ia(i) ),'_J', ...
                                      num2str( ja(i) ) ) ) ); 
        end
        
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
                
                if pltflag_regression == true
                    figure             
                    PLR = plotregression( dataDRT{j,1}, dataDRT{j,2}, ...
                        strcat('(',num2str( ia(i) ),    ...
                                   ',',num2str( ja(i) ),')',...
                                ' DRT=',num2str( drt(j) ),' tol_s=',num2str(seps) ) );
                    
                    setLRPlot(PLR); % graph appearance
                end
                
                % csv file [ logphiz, logRQI, depth ] per DRT                
                fprintf('Exporting CSV file of logs/depth data for well(%d,%d); DRT %d... \n', ia(i), ja(i), drt(j) );  
                fname = fullfile( '../csv/', strcat(regname,num2str( ia(i) ),'_J', num2str( ja(i) ),'_DRT_',num2str( drt(j) ),'_LogsDepth','.csv' ) );
                auxmat = [ dataDRT{j,1} dataDRT{j,2} dataDRT{j,3} ];
                                
                % header
                head = {'LogPhiz,';'LogRQI,';'z'}; 
                head = head';
                txt = sprintf('%s\t',head{:}); 
                txt(end) = '';                
                dlmwrite(fname,txt,'');                                
                dlmwrite(fname,auxmat,'-append');                                

                % csv file: fit data
                fprintf('Exporting CSV file of regression fit data for well(%d,%d); DRT %d... \n', ia(i), ja(i), drt(j) );  
                fname = fullfile( '../csv/', strcat(regname,num2str( ia(i) ),'_J', num2str( ja(i) ),'_DRT_',num2str( drt(j) ),'_FitData','.csv' ) );
                
                % REMARK: rj value from regression() above is computed by Matlab 
                %         as it is. To have R-squared (R^2) value 
                %         (Pearson correlation), we should save rj*rj, i.e.
                %         ----> dataDRT{j,4}*dataDRT{j,4}
                %
                aux = [ dataDRT{j,4}*dataDRT{j,4} dataDRT{j,5} dataDRT{j,6} ];
                
                % header
                head = {'R^2,';'slope,';'offset'}; 
                head = head';
                txt = sprintf('%s\t',head{:}); 
                txt(end) = '';                
                dlmwrite(fname,txt,'');                                
                dlmwrite(fname,aux,'-append');                                
                
                % print to file   
                if pflag_regression == true
                    print('-dpdf','-r0',fullfile( '../figs/', ...
                       strcat(regname,num2str( ia(i) ),'_J', ...
                                      num2str( ja(i) ), ...
                                      '_DRT_',num2str( drt(j) ) ) ) );                
                end
                
            end

        end   
                                
        if pltflag_HFULoc == true
            % Selecting the best DRTs for plotting
            depths = cell(1,length(drtgood));
            for k = 1:length(drtgood)
                id = drtgood(k);
                idd = find( drt(:) == id );
                depths{1,k} = dataDRT{ idd,3 };
            end

            % plot the best HFU locations per well
            %drt2color = 13;
            %plotHFULoc( drtgood,depths,ia(i),ja(i),drt2color ); 
        end
                

    end

end

%% ENDING
d.printings(d.progStat{2});
dm.deactivateLog;

