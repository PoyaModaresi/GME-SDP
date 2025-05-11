function MUBS = GenerateMUB(n,d)
% This function generates all possible MUBs
% - n: number of parties
% - d: dimension of local subsystems

%
B=mub(d);
amount=size(B,3);

%Generate every combo of k with repeated ind
combo=GenAllK(n,d);

%Make the d^n X amount matrix where amount is the amount of MUBS available
%and d^n is the amount of projectors/mesurment per MUBS.
shell=cell(d^n,amount);



%First loop all the mub
for i = 1:amount
    basisMatrix=B(:,:,i);
    for l =1:size(combo,1)
        %Loop the vector for the current MUB
        currentK=combo(l,:); %The current vector
        vector=1;
        for j =1:length(currentK)
            vector = Tensor(vector,basisMatrix(:,currentK(j))); %Tensor them into a state
        end
        shell{l,i}= vector*vector'; %make the proj mesurment

    end
end

MUBS=shell;


