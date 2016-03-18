classdef SPEDirManager < handle
    %SPEDirManager directory management    
    
    properties
    end
    
    methods (Static)
        
        % ACTIVATELOG
        function activateLog(mainfile)  
            if exist('../log','dir') ~= 7; mkdir('../log'); end
            delete( strcat('../log/',mainfile,'.log') ); % delete old file
            diary( strcat('../log/',mainfile,'.log') ); diary on
        end
        
        % DEACTIVATELOG
        function deactivateLog
            diary off
        end
                                
    end
    
end

