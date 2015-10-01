function execCSVCondenser
% !!! Under development
    
disp('Executing CSV condenser...');

    csvdir = '../csv/';

    outdir = '../csv/condensed';

    fd = dir(fullfile(csvdir,'*LogsDepth.csv'));

    if ~isempty(fd)

        nw = length(fd); % number of files 

        Aout = [];
        for i = 1:nw
            fdname = fd(i).name;        

            [ini,eni] = regexp(fdname,'I\d*_');      % find I index in filename
            wI = str2double( fdname(ini+1:eni-1) ); 

            [inj,enj] = regexp(fdname,'J\d*_');      % find J index in filename
            wJ = str2double( fdname(inj+1:enj-1) ); 

            [int,ent] = regexp(fdname,'DRT_\d*_');   % find DRT in filename
            wDRT = str2double( fdname(int+4:ent-1) ); 

            fprintf('Condensing CSV file of DRT=%d for well (%d,%d)...\n',wDRT,wI,wJ);
            %{
                 %TODO Improve the way to insert a separation header string

                Writing CSV file ( Mandrake example ) 
                =====================================

                  1   1   1  } Header for, let' say, DRT = 1
                  1   1   1  }        

                [.   .   .]  } contents for DRT = 1          

                  4   4   4  } Header for, let' say, DRT = 4
                  4   4   4  }        

                [.   .   .]  } contents for DRT = 1     , 

                                        where [. . .] = [ logphi logDRT k ]


            %}
            A = csvread( fullfile( csvdir, fdname ) );         
            A = [ wDRT*ones(2,3); A ];

            Aout = [ Aout; A ];

            outname = strcat( 'BestDRTsCondensed_I',num2str(wI),'_J_',num2str(wJ) );         
            csvwrite( strcat( fullfile(outdir,outname),'.csv'),Aout);

        end

    else
        error('csv dir is empty.');
    end

end


