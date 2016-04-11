function depaths = addDependencies( varargin )
%ADDDEPENDENCIES add path to dependent files

for i=1:length(varargin)
    depaths{i} = varargin{i};
end

