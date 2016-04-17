function cleanAllSPE
%CLEANALLSPE delete post-processed files
    
    dnames = {'../csv/' ,...
              '../figs/',...
              '../log/' ,...
              '../vtk/' ,...
              '../tmp/' ,...
              '../txt/' ,...
              '../img/'      };     % post-processing dirs
    
    for dn = 1:length(dnames)                
        
        delete( strcat(dnames{dn},'*.*') ); % common files        
        
        dr = dir(dnames{dn});   % remove subdirs recursively
        if length(dr) > 2
            for j = 3:length(dr)            
                fname = fullfile(dnames{dn},dr(j).name);
                if exist(fname,'file') == 7 
                    rmdir(fname,'s' );
                else
                    delete(fname);
                end
            end
        end
        
    end
           
    
    
    
end

