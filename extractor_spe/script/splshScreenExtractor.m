function splshScreenExtractor( )

disp( repmat('-', [1 49] ) );
disp('           SPE 2 - RESERVOIR EXTRACTOR           ');
disp( repmat('-', [1 49] ) );
disp(' ');
disp('By:    Dr. Gustavo Peixoto de Oliveira');         
disp('       Dr. Waldir Leite Roque');
disp(' ');
disp('       @Federal University of Paraiba - Brazil');
disp(' ');
disp( repmat('*', [1 54] ) );
fprintf(['ATTENTION! This script is based on the synthetic   \n',...
      'field SPE Model 2 whose dataset is available online.  \n',... 
      'To work, the path to the reservoir input .dat files   \n',... 
      'of porosity and permeability must be provided. If     \n',...
      'you are using it at first time and do not have the    \n',... 
      '.mat files already computed, choose [no] in the rerun \n',... 
      'prompt. Otherwise, choose [yes] to bypass assembling  \n',... 
      'operations and additional time-consuming computations.\n']);       
disp( repmat('*', [1 54] ) );
disp('S T A R T I N G     P R O G R A M');
disp(' ');

end

