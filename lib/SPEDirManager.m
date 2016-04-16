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
        
        function [fout,wname] = createWellDir(ic,jc,basedir)
            % CREATEWELLDIR creates dir for specified well
            %   input: 
            %           ic,jc: well's surface coordinates
            %         basedir: dir where to create ('mat','txt',etc.)
            %  output:
            %            fout: local path to output dir
            %           wname: well name
            %
            wname = strcat('Well_I',num2str(ic),'_J',num2str(jc));                        
            frel = strcat('../',basedir);
            fout = strcat(frel,'/',wname);            
            if exist(frel,'dir') ~= 7; mkdir(frel); end    
            if exist(fout,'dir') ~= 7; mkdir(fout); end    
            fout = strcat(fout,'/');            
        end
        
        function vfile = setVOIFile(fout,wname,drt)
            % SETVOIDIR set VOI file per DRT and well
            %
            %   input: 
            %           fout: output from CREATEWELLDIR;
            %          wname: well name
            %            drt: drt value (int)
            %   output:
            %           vfile: VOI file
            vfile = fullfile(fout,...
                         strcat('VOI_DRT_',num2str(drt),'_',wname,'.mat'));
            
        end
        
        function zipCmd(scal,flag)    
            % ZIPCMD calls zip to compress images 
            %
            %   input:
            %
            %           scal: scalar name('phi','kx,'ky,'kz','drt');
            %           flag: execute zip (true,false)            
            out = strcat('../img/',scal);
            in = strcat('../img/',scal,'_img/');
            cmd = sprintf('zip -r %s %s',out,in);            
            if flag == true; 
                [~,~] = system(cmd,'-echo'); 
            end            
        end
        
        function createDirStructure
            
            if exist('../csv/','dir') ~= 7; 
                mkdir('../csv/');                 
            end    
            
            if exist('../vtk/','dir') ~= 7; 
                mkdir('../vtk/');                 
            end    
            
            if exist('../txt/','dir') ~= 7; 
                mkdir('../txt/');                 
            end    
            
            if exist('../mat/','dir') ~= 7; 
                mkdir('../mat/');                 
            end    
            
            if exist('../img/','dir') ~= 7; 
                mkdir('../img/');                 
            end    
            
            if exist('../figs/','dir') ~= 7; 
                mkdir('../figs/');                 
            end    
            
            if exist('../tmp/','dir') ~= 7; 
                mkdir('../tmp/');                 
            end    
            
            if exist('../log/','dir') ~= 7; 
                mkdir('../log/');                 
            end    
            
        end
        
    end
    
end

