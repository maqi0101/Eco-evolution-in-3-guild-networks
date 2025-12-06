% This function is used to calculate functional diversity (Rao's quadratic entropy)

function FDQ = get_FDQ(trait,ymea)
FDQ=0;
S=length(ymea);   % total number of species
p=ymea/sum(ymea); % relative abundance
for i=1:S-1
    for j=i+1:S
        FDQ=FDQ+abs(trait(i)-trait(j))*p(i)*p(j);
    end
end

