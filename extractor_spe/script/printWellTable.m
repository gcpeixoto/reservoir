function printWellTable( ia, ja, N )
%{
    PRINTWELLTABLE print table with well's surface coordinates.

%}
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