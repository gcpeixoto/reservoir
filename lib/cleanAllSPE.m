function cleanAllSPE
    delete( strcat('../csv/*.csv') );
    delete( strcat('../csv/subdomain/*.csv') );
    delete( strcat('../log/*.log') );
    delete( strcat('../vtk/*.vtk') ); 
    delete( strcat('../tmp/*.*') );     
    delete( strcat('../figs/subdomain/*.*') );
    delete( strcat('../figs/graphPath/*.*') );
    delete( strcat('../dat/subdomain/*.dat') );
    
end

