function L = reductionMap(rho,k,n,d,cut,perm)
% This function applies the reduction map
% - rho: state density matrix
% - k: GME-dimension
% - n: number of parties
% - d: dimension of local subsystems
% - cut: set on which apply the k-positive operator \Lambda
% - perm: permutation vector to restore the order of the parties after the
% partial trace


%Identity
id = eye(d);
a = length(cut);
idA = eye(d^a);
idn = Tensor(id,n);

%Generalised reduction map
L = PermuteSystems(Tensor(idA,PartialTrace(rho,cut,d*ones(1,n))),perm) - 1/k*rho;

end