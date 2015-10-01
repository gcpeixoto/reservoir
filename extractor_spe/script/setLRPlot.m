function setLRPlot( hdl )
%{
    SETLRPLOT Modifies plotregression appearance
    
    The plotregression function returns the handle to a figure. 
    This figure has 3 children: legend, axis, and a uicontrol. 
    For a simple call the uicontrol is not visible. 
    The axis has also has 3 children: data, fit, y = T. 
    To get what you want we need to delete the third child 
    of the second child and change the marker of the 
    first child of the second child. We then need to regenerate 
    the legend since it does not dynamically update.
    
    From: http://stackoverflow.com/questions/18770166/matlab-plot-regression-function
%}
h = get(hdl, 'Children');
hh = get(h(2), 'Children');
delete(hh(3))
set(hh(1), 'Marker', 'o','MarkerFaceColor','k','MarkerEdgeColor','k');
legend('FZI Fit','Points','Location', 'NorthEastOutside');


end

