function drtVec = getDRTFromFileName( matFiles )
%GETDRTFROMFILENAME extracts DRT values from file names
%            input: .mat files (struct)
%           output: DRT values (nx1)

% test
assert(isstruct(matFiles),'getDRTFromFileName > argument is not a struct');

% gets DRT values
nf = length(matFiles);
drtVec = zeros(nf,1);
for i = 1:nf
    [~,aux,~] = fileparts(matFiles(i).name);    
    val = regexp(aux,'(_)\d*(_)','match');
    val = val{1}(2:end-1); % only number
    val = str2num(val); % converts
    drtVec(i) = val;    
end


