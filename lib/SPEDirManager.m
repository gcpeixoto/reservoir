classdef SPEDirManager < handle
    %SPEDirManager directory management    
    
    properties
    end
    
    methods (Static)
                
        function activateLog(mainfile)  
            % ACTIVATELOG open diary to create a log file for the running
            %             file
            
            if exist('../log','dir') ~= 7; mkdir('../log'); end
            delete( strcat('../log/',mainfile,'.log') ); % delete old file
            diary( strcat('../log/',mainfile,'.log') ); diary on
        end
                
        function deactivateLog
            % DEACTIVATELOG close log file
            diary off
        end
        
        function fout = createWellDir(ic,jc,basedir)
            % CREATEWELLDIR creates dir for specified well
            %   input: 
            %           ic,jc = well's surface coordinates
            %         basedir = dir where to create ('mat','txt',etc.)
            %  output:
            %            fout = local path to output dir
            %
            wname = strcat('Well_I',num2str(ic),'_J',num2str(jc));                        
            frel = strcat('../',basedir);
            fout = strcat(frel,'/',wname);            
            if exist(frel,'dir') ~= 7; mkdir(frel); end    
            if exist(fout,'dir') ~= 7; mkdir(fout); end    
            fout = strcat(fout,'/');
            
        end
                                
    end
    
end

