classdef SPEDisplay < handle
    %SPEDISPLAY class responsible for visual and display formatting
    
    properties 
                
        author1 = 'Dr. Gustavo Peixoto de Oliveira';
        author2 = 'Dr. Waldir Leite Roque';
        inst = 'Federal University of Paraiba';
            
        % execution status
        progStat = {'S T A R T I N G  . . .',...
                    'F I N I S H I N G  . . .'};        
        % messages 
        
        % mainExtractorSPE
        extSPERerun   = '----> Rerun extractor? [0] no; [1] yes \n';
        extSPERebuild = '----> Rebuilding porosity and permeability matrices...';
        extSPEReload  = '----> Reloading porosity and permeability matrices...';
        extSPEWExt    = ['----> Choose method for well extraction:',... 
                         '[0] N random wells; [1] specific well: \n'];                     
        extSPEChooseN = '----> Choose number of random wells to extract: \n';
        extSPEWellI   = '-------> Enter with I coordinate of the well: \n'
        extSPEWellJ   = '-------> Enter with J coordinate of the well: \n'
        extSPEGrossCSV= '----> Export gross well data to CSV? [0] no; [1] yes. \n';
        extSPELinRegr = '----> Compute histograms and regression analysis? [0] no; [1] yes. \n';
    
        % mainExtractorSPESubvolume
        extSPESP      = '-----> Choose P ring radius: \n';
        
        % mainVolume2Image
        vol2Im = '-----> Compressing images to .zip...';
        
    end
    
    methods (Static = true) % not-dependent on class objects
        %% PRINTS 
                
        function printings(varargin)                        
         % PRINTINGS sequential on-screen prints
            for i = 1:length(varargin); disp(varargin{i}); end
        end
            
        
        function printSplScreen(mainfile)                        
        % PRINTSPLSCREEN execution filename print
            tit = strcat( 'SPE EXTRACTOR::',upper(mainfile) ); 
            nchar = length(tit);
            hline = repmat('-',[1, round(0.2*nchar)]);
            hline = strcat(hline,repmat('-',[1,nchar]),hline,'\n');
            htit = repmat(' ',[1, round(0.2*nchar)]);            
            fprintf(hline);
            fprintf('%s %s %s \n',htit,tit,htit);            
            fprintf(hline);                                     
        end
        
        function printWellTable( ia, ja, N )
        % PRINTWELLTABLE print table with well's surface coordinates.

            disp('');
            disp('---- Well Table ----');
            fprintf('Well\t I\t J\t\n'); % header

            i = 1;
            while i <= N % wells
                fprintf('%d\t %d\t %d\t\n', i, ia(i), ja(i) );    
                i = i + 1;
            end
            disp( repmat('-', [1 20]) );

        end
        
        % TODO mainExtractorSPE
        function printWSatTable( min1,max1,min2,max2,min3,max3 )
        % PRINTWSATTABLE print table with well's surface coordinates.
            disp('');
            disp('---- Irreducible water saturation bounds ----');
            fprintf('Layer\t wsiA\t wsiB \n'); % header
            fprintf('1\t %g\t %g \n', min1, max1 );    
            fprintf('2\t %g\t %g \n', min2, max2 );    
            fprintf('3\t %g\t %g \n', min3, max3 );    
            disp( repmat('-', [1 45] ) );

        end
                      
        %%  WARNINGS
                
        function extractorSPEDependency        
            warning('this method requires .mat files from mainExtractorSPE.m.');
        end
              
        
        function graphDataDependency            
            warning('this method requires .mat files from mainDRTGraphData.m');
        end
                
        
        function VOIgraphDataDependency
            warning('this method requires .mat files from mainVOIDRTGraphData.m');
        end
           
        
        function VOIgraphMetricsDependency
            warning('this method requires .mat files from mainVOIDRTGraphMetrics.m');
        end
        
        
        function VOIConnectionsDependency
            warning('this method requires .mat files from mainVOIConnections.m');
        end
        
        
        function dispTest
            warning('THIS METHOD IS UNDER TESTING!');
        end    
                
        
        function dispMsg(msg)
            disp(msg);
        end
                        
        function dispNotValid
            error('Option not valid!');
        end
        
        % DISPCCORD
        % mainExtractorSPESubvolume
        function msg = dispCCoord(coord)
                        
            switch coord
                case 1
                    c = 'I';
                case 2
                    c = 'J';
                case 3
                    c = 'K';                    
            end
            msg = strcat('-----> Choose central voxel coordinate::',c,'\n');
        end                
        
        function setOptions                                    
            alw = 1.0; % AxesLineWidth    
            afs = 14;  % AxesFontSize
            set(0,'DefaultAxesFontSize',afs);          
            set(0,'DefaultAxesLineWidth',alw);
            set(0,'DefaultAxesFontUnits','points');
            set(0,'DefaultAxesGridLineStyle',':');
            set(0,'DefaultAxesFontName','Helvetica');
            set(0,'DefaultTextFontName','Helvetica');
            set(0,'DefaultAxesTickDir','out');
            format long
        end                 
            
    end
    
end

