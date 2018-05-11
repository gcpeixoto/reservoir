function plotHFULoc( DRT, Kd, I, J )
% PLOTHFULOC Plot HFU locations per well based on the best DRTs

figure
xlabel('$ DRT $','interpreter','latex');  
ylabel('$ z $','interpreter','latex');
title( strcat('High-performance DRTs - (',num2str(I),',',num2str(J),')' ) );            
grid on, box on
set(gca,'YDir','reverse');

for i = 1:length(DRT)    
    hold on    
    xlim( [ min(DRT-1) max(DRT+1) ] );    
    plot(DRT(i), fliplr( Kd{i} ),'ok','MarkerFaceColor','k');         
end

end

