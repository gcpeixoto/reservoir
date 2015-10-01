function saveExpurgatedIds( idvec, i, j )
%{
    SAVEEXPURGATEDIDS Write to file the reservoir's 
     depth z-coordinates where a value of 0 was
     found for porosity. These values are removed in the analysis
     and this file is only intended to identify such points.
%}

for k=1:length(idvec)    
    aux = [ i j idvec(k) ];
    fname = fullfile( '../csv/',strcat( 'ExpurgatedPoints_I_',num2str(i),'J_',num2str(j),'.csv' ) );    
    csvwrite(fname, aux);
end

