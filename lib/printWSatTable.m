function printWSatTable( min1,max1,min2,max2,min3,max3 )

disp('');
disp('---- Irreducible water saturation bounds ----');
fprintf('Layer\t wsiA\t wsiB \n'); % header
fprintf('1\t %g\t %g \n', min1, max1 );    
fprintf('2\t %g\t %g \n', min2, max2 );    
fprintf('3\t %g\t %g \n', min3, max3 );    
disp( repmat('-', [1 45] ) );

end