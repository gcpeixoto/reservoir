function cleanAllSPE
%CLEANALLSPE delete post-processed files

    delete( strcat('../csv/*.csv') );    
    delete( strcat('../log/*.log') );
    delete( strcat('../vtk/*.vtk') ); 
    delete( strcat('../tmp/*.*') );         
    delete( strcat('../figs/*.*') );
    delete( strcat('../txt/*.*') );
    delete( strcat('../dat/subdomain/*.dat') );
    
end

