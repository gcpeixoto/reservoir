classdef SPEDisplay
    %SPEDISPLAY class responsible for visual and display formatting
    
    properties 
        
        % execution status
        progStat = {'S T A R T I N G  . . .',...
                    'F I N I S H I N G  . . .'};        
        author1 = 'Dr. Gustavo Peixoto de Oliveira';
        author2 = 'Dr. Waldir Leite Roque';
        inst = 'Federal University of Paraiba';
        
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
    
    methods (Static)
        
        % PRINTINGS sequential on-screen prints
        function printings(varargin)                        
            for i = 1:length(varargin); disp(varargin{i}); end
        end
            
        % PRINTSPLSCREEN execution filename print
        function printSplScreen(mainfile)                        
            tit = strcat( 'SPE EXTRACTOR::',upper(mainfile) ); 
            nchar = length(tit);
            hline = repmat('-',[1, round(0.2*nchar)]);
            hline = strcat(hline,repmat('-',[1,nchar]),hline,'\n');
            htit = repmat(' ',[1, round(0.2*nchar)]);            
            fprintf(hline);
            fprintf('%s %s %s \n',htit,tit,htit);            
            fprintf(hline);                                     
        end
                        
        % EXTRACTORSPEWARNING
        function extractorSPEwarning
            warning(['running mainExtractorSPE for the first time? ',...
                     'Choose rerun (option [1]) below to load .mat files.']);            
        end
        
        % EXTRACTORSPEDEPENDENCY
        function extractorSPEDependency
            warning('this method requires .mat files from mainExtractorSPE.m.');
        end
        
        % GRAPHDATAPEDEPENDENCY
        function graphDataDependency
            warning('this method requires .mat files from mainDRTGraphData.m');
        end
        
        % VOIGRAPHDATAPEDEPENDENCY
        function VOIgraphDataDependency
            warning('this method requires .mat files from mainVOIDRTGraphData.m');
        end
        
        % GRAPHMETRICSPEDEPENDENCY
        function VOIgraphMetricsDependency
            warning('this method requires .mat files from mainVOIDRTGraphMetrics.m');
        end
        
        % DISPMSG
        function dispMsg(msg)
            disp(msg);
        end
                 
        % DISPNOTVALID
        function dispNotValid
            error('Option not valid!');
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
        
        function dispTest
            warning('THIS METHOD IS UNDER TESTING!');
        end     
            
    end
    
end

