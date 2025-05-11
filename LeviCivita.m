function e=LeviCivita(L)
% This function calculates the Levi-Civita function
% - L: Permutation

UL=unique(L);

if isempty(L)==1
    warning("Please input non empty list")
elseif length(L)-length(UL)>0
    e=0;
else
    S=1;
     for i = 1:length(L)
         for j = i:length(L)
             k=L(j)-L(i);
             if k~=0
                 S=S*sign(k);
             end
         end
     end
     e=S;
end





