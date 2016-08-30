function [ ncolfile ] = saveNcolFile( Madj, w, varargin )
% SAVELGLEFILE write adjacency matrix edge list to a .NCOL file 
%              to interface with Python iGraph.
%  
% description:
%   NCOL format is a simple edge list determing the connectivity 
%   of a graph (See http://lgl.sourceforge.net/#FileFormat)
%   
%   input: 
%              Madj: adjacency matrix (sparse nxn)
%                 w: weights of the edges, if any (nx1)

% temporary edge file
ncolfile = '../tmp/edges.ncol';

    
if issparse(Madj)
    iline = [];        
    [i,j] = find(Madj == 1);
    edges = [i,j];
    
    for i = 1:size(edges,1)
        for j = i+1:size(edges,1)
            if (i ~=j && all(edges(i,:) == fliplr(edges(j,:))))
                iline = [iline; i];
            end
        end
    end
    
    if nargin == 1
        edges = edges(iline,:);
    else
        edges = [edges(iline,:),w(iline)];        
    end
    
    dlmwrite(ncolfile,edges,'delimiter',' '); % writing file

else
    error('Matrix is not sparse');
end
    