function maxflow = getFlowData(capfile,ssfile)
% GETMETRICSFLOW get flow data computed from NETWORKX.

pyi = setPyInterpreter('/usr/local/bin/python');
pycmd = strcat(pyi, ' ../py/flow.py');        
[status, maxflow] = system(pycmd);

maxflow = str2double( maxflow(17:end) );

if status ~= 0
    warning('Check flow data calculation.')
end

% delete temporary files
%delete(capfile); 
%delete(ssfile); 

end