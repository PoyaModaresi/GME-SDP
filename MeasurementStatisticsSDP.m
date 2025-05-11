function v = MeasurementStatisticsSDP(n,d,k,state,option,channels,MUBS)



    fprintf('Initializing SDP problem for n=%d, d=%d, k=%d...\n', n, d, k);
%-------------------------------------------------------------------------%
%This function evaluates the measurment statistic SemidefiniteProgram (SDP)
%for a target n-partite d-dimensional state using YALMIP
%Inputs:
% - n: number of parties
% - d: dimension of local subsystems
% - k: GME-dimension
% - state: pure target density matrix
% - option: application mode of the generalised reduction map
% -- 0: Partial trace on the smallest subset
% -- 2: Partial trace on the biggest subset
% -- 1: Map applied on both subsets
% - channels: Apllication mode for the quanum channel
% -- 0 Depolarizing Channels
% -- 1 Dephasing  Channels
% -- MUBS: Number of MUBs. If not given it will do all possible MUBS
% available

%Output:
% - v: maximum visibility at which a noisy state has GME-dimension = k
%-------------------------------------------------------------------------%

%pure density state
rho=state*state';

%identity
id=eye(d);

%amount of possible bipartation
s=2^(n-1)-1;

%list of possible bipartation
bipartitions = SetPartition(n,2);

%initilize constraints list and visability variable
constraints=[];
v = sdpvar(1);


P = GenerateMUB(n, d);
num_proj = size(P, 1);
MaxMub=size(P, 2);
if nargin<7
    num_mubs = MaxMub;
else
    if MUBS<=MaxMub
        num_mubs = MUBS;
    else
         warning('To many MUBS requested, please pick' + MaxMub + "or less")
    end
end

%initilize bipartation states and sum of all bipartation states
    fprintf('Creating SDP variables...\n');
SigmaTotal=0;
for i = 1 : s
    Sigma{i} = sdpvar(d^n,d^n,'hermitian','complex');
    %    SigmaTotal = SigmaTotal + Sigma{i};
    constraints=[constraints, Sigma{i}>= 0,];
end
    SigmaTotal = sum(cat(3, Sigma{:}), 3);

%Apply reduction map
    fprintf('Applying reduction map...\n');
for i = 1 : s
    %Discriminate between smallest and biggest subsets
    b = cell2mat(bipartitions{i}(1));
    c = cell2mat(bipartitions{i}(2));
    if length(b) <= length(c)
        set = [b];
        comp = [c];
    else
        set = [c];
        comp = [b];
    end
    %Map on the smallest subset
    if option == 0
        O="Smallest Subset";
        order = [set comp];
        [~,perm] = sort(order);
        reduceMap = reductionMap(Sigma{i},k,n,d,set,perm);
        constraints = [constraints, reduceMap >= 0];
    %Map on the both subsets
    elseif option == 1
        O="Both Subsets";
        order1 = [set comp];
        [~,perm1] = sort(order1);
        order2 = [comp set];
        [~,perm2] = sort(order2);
        reduceMap1 = reductionMap(Sigma{i},k,n,d,set,perm1);
        reduceMap2 = reductionMap(Sigma{i},k,n,d,comp,perm2);
        constraints = [constraints, reduceMap1 >= 0, reduceMap2 >= 0];
    elseif option == 2
    %Map on the bigger subset
        O="Bigger Subset";
        order = [comp set];
        [~,perm] = sort(order);
        reduceMap = reductionMap(Sigma{i},k,n,d,set,perm);
        constraints = [constraints, reduceMap >= 0];
    else
        error('option not valid: Please pick between 0, 1 and 2');
    end
end



%Constraints for the choosen quantum channel
    fprintf('Applying quantum channel constraints...\n');
if channels ==0
    Ctext="Depolarizing";
    inside=v*rho + (1-v)*Tensor(id,n)/d^n;
    for j = 1:num_mubs
        for i = 1:num_proj
        M = P{i, j};
        constraints = [constraints,trace(SigmaTotal * M) == trace((inside) * M)];
        end
    end
elseif channels==1
        Ctext="Dephasing";
    comp = cell(1, d);
    for i = 1 : d
        comp{i} = id(:,i);
    end

    Dephase=0;

    for i = 1:d
        iState=Tensor(comp{i}, n);
        Dephase=Dephase+iState*iState';
    end
    inside=v*rho + (1-v)/d*Dephase;

    for j = 1:num_mubs
        for i = 1:num_proj
        M = P{i, j};
        constraints = [constraints,trace(SigmaTotal * M) == trace((inside) * M)];
        end
    end
else
    error("Channel not valid: Please pick between 0 and 1")
end




%Solve using YALMIP + MOSEK
    fprintf('Solving SDP...\n');
ops = sdpsettings('solver', 'mosek', 'verbose', 1, 'cachesolvers', 1, ...
                  'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3, ...  % Increase tolerance
                  'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-3, ...
                  'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 600);  % 10-minute timeout

diagnostic=solvesdp(constraints,-v,ops);

    if diagnostic.problem ~= 0
        warning('Solver failed with status: %d', diagnostic.problem);
        v = NaN;
        return;
    end

v=double(v);

%display result
disp("For " + Ctext + " channel the critical visibility when applying the positive map to " + O + " is " + v)
end