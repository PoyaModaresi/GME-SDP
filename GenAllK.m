function  combos= GenAllK(n,d)
% This function generates all iterations of n numbers ranging from d
% - n: number of parties
% - d: dimension of local subsystems

v=1:d;
indexing=cell(1,n);
[indexing{:}] = ndgrid(v);
combos = cell2mat(cellfun(@(x) x(:), indexing, 'UniformOutput', false));
