function plotHFULoc( DRT, Kd, I, J )
%plotHFULoc Plot HFU locations per well based on the best DRTs

figure
xlabel(' $ DRT $','interpreter','latex');  
ylabel('depth');
title( strcat('HFU locations - (',num2str(I),',',num2str(J),')' ) );            
grid minor
set(gca,'YDir','reverse');
fname = 'HFULocationsI';

for i = 1:length(DRT)    
    hold on    
    xlim( [ min(DRT-1) max(DRT+1) ] );    
    plot( DRT(i), fliplr( Kd{i} ),'ko', 'MarkerFaceColor', 'k' );         
end

% print to file    
print('-dpdf','-r0',fullfile( '../figs/', ...
          strcat(fname,num2str( I ),'_J', num2str( J ) ) ) );                

end

