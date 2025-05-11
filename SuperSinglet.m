function state = SuperSinglet(d)
% This function creates the supersinlget state
% - d: dimension of the supersinglet state


%List of numbers from 0 to d-1
number= 0:d-1;

%All possible permutations of number
C = perms(number);

%Normalazation factor
norm=1/sqrt(factorial(d));

%Identity
id = eye(d);

%Computational basis
comp = cell(1, d);
for i = 1 : d
    comp{i} = id(:,i);
end

%Making supersinglet state
state=0;
for i = 1 : length(C)
    perm=C(i,:);
    e=LeviCivita(perm);
    tensor_input = cell(1, d);
    for j = 1:d
        tensor_input{j} = comp{perm(j) + 1};
    end
    

    state = state + e * norm * Tensor(tensor_input);
    %state=sparse(state);
    
end


%Quick test
%012-021-102+120+201-210
%123-132-213+231+312-321
%test=norm*( Tensor(comp{1},comp{2},comp{3}) - Tensor(comp{1},comp{3},comp{2}) - Tensor(comp{2},comp{1},comp{3}) +Tensor(comp{2},comp{3},comp{1}) + Tensor(comp{3},comp{1},comp{2})-Tensor(comp{3},comp{2},comp{1}) );

end