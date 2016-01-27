function [ outfile ] = saveCapacityTable( Madj, c )
% SAVECAPACITYTABLE write adjacency matrix edge list to file 
%   
%   input: 
%           Madj: adjacency matrix
%              c: capacity values

% temporary edge file
outfile = '../tmp/capacity';

if issparse(Madj)
    [i,j] = find(Madj == 1);
    deltac = c(i) - c(j); % DELTA: pressure gradient
    capacity = [i,j,deltac];
    dlmwrite(outfile,capacity,'delimiter',' '); % writing file
else
    error('Matrix is not sparse');
end