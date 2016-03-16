classdef SPEDirManager
    %SPEDirManager directory management    
    
    properties
    end
    
    methods (Static)
        
        % ACTIVATELOG
        function activateLog(mainfile)    
            delete( strcat('../log/',mainfile,'.log') ); % delete old file
            diary( strcat('../log/',mainfile,'.log') ); diary on
        end
        
        % DEACTIVATELOG
        function deactivateLog
            diary off
        end
                                
    end
    
end

