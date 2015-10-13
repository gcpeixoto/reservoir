function [PMax,IC,JC,KC] = getMaxRing( I,J,K )
%{ 
    getMaxRing - gets the maximum ring value 
                 allowed for the grid
    
    input: 
        I,J,K: grid bounds        

    output:      
         PMax: maximum ring value allowed for the grid
     IC,JC,KC: central seed coords
        
%}

%{
    seed location
    If the grid has some even bound, there's no 
    an integer median. Otherwise, the median
    can be used. 
%}
if mod(I,2) == 0 || mod(J,2) == 0 || mod(K,2) == 0 % even 
    IC = floor( median(1:I) );
    JC = floor( median(1:J) );
    KC = floor( median(1:K) );

    % maximum ring value allowed for the grid
    PMax = min( [I-IC J-JC K-KC] ) - 1;
else % odd
    IC = median(1:I);
    JC = median(1:J);
    KC = median(1:K);

    % maximum ring value allowed for the grid
    PMax = min( [I-IC J-JC K-KC] );
end


end

