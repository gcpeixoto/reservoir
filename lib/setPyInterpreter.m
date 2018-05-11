function str = setPyInterpreter( path )
% SETPYINTERPRETER sets the Python interpreter to use
%
% input: 
%   path: full path to interpreter
%
%
% REMARK: note that call `which python` from bash may differs 
%         from `system('which python') called from inside MATLAB.
%         If so, the parsing of python scripts that import non-standard
%         python modules will not work.
%
% See discussions here: 
%
%   <http://goo.gl/SWQtcx> and <http://goo.gl/4oQdlM>
%

if isempty(path) || ~ischar(path)
    % TODO
    % Add if-clause to included a pydefault for Windows, Linux
    pydefault = '/usr/local/bin/python'; % MAC OS
    warning( strcat('Setting python interpreter to: ',pydefault) );
    str = pydefault;
else
    str = path;
end

