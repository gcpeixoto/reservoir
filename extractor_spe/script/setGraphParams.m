function [ capfile, ssfile ] = setGraphParams( Madj, c, src, sink )
% SAVEGRAPHPARAMS write info to files in order to create interface
%                 for graph computation
%   
%   Parameters exported to file: edge list of adjacency matrix
%                                edge weights
%                                source and sink node
%   
%   input: 
%           Madj: adjacency matrix
%              c: capacity values (weights for each edge)
%            src: source node for maximum flow problem
%           sink: sink node for maximum flow problem
%
%   output:
%       capfile: file containing source | target | capacity
%        ssfile: file containing sink, source node
%


% temporary edge file with weigths
capfile = '../tmp/capacity';

% temporary source,sink file
% source node | sink node
ssfile = '../tmp/sinksource';

if issparse(Madj)
    [i,j] = find(Madj == 1);
    deltac = c(i) - c(j); % DELTA: pressure gradient
    capacity = [i,j,deltac];            
    dlmwrite(capfile,capacity,'delimiter',' '); % writing file
else
    error('Matrix is not sparse');
       
end


dlmwrite(ssfile,[src,sink],'delimiter',' '); % writing to file

end