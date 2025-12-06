% This function is used to update interspecies interactions
% (Update the interspecific trait matching, and rewiring one link)

function [match_PH,match_PM]=updata_interaction(Sp,Sh,Sm,trait,index_P2A,match_PM,match_PH,sigma,eta)

% species traits after coevolution
Z_P=trait(1:Sp);             % plant
Z_H=trait(Sp+1:Sp+Sh);       % herbivore
Z_M=trait(Sp+Sh+1:Sp+Sh+Sm); % pollinator

%%
% update the interspecific trait matching
for i=1:Sp
    for j=1:Sh
        if match_PH(i,j)>0
            match_PH(i,j)=trait_match(Z_P(i),Z_H(j),sigma);
        end
    end
end

for i=1:Sp
    for j=1:Sm
        if match_PM(i,j)>0
            match_PM(i,j)=trait_match(Z_P(i),Z_M(j),sigma);
        end
    end
end

%%
% calculate the difference of trait between animals and their potential plant partners
TD_HP=zeros(Sh,Sp);
for i=1:Sh
    for j=1:Sp
        if match_PH(j,i)==0  % Judge whether plant j is a potential partner of herbivore i
            TD_HP(i,j)=abs(Z_H(i)-Z_P(j)); % difference of trait
        end
    end
end

TD_MP=zeros(Sm,Sp);
for i=1:Sm
    for j=1:Sp
        if match_PM(j,i)==0
            TD_MP(i,j)=abs(Z_M(i)-Z_P(j));
        end
    end
end

%%
% rewiring one link

Type=randi(2);   % rewiring in H or M guild  (1:H->P; 2:M->P)

if (Type==1)  % rewiring H->P: herbivore to plant

    idH=randi(Sh); % randomly select a herbivore
    while isempty(find(match_PH(:,idH),1))  % ensure the herbivore has at least one link
        idH=randi(Sh);
    end
    
    idP_old=index_P2A(idH);  % the ID number of the plant (P_j) with the least evolutionary effect on H_i
    TD_old=abs(Z_H(idH)-Z_P(idP_old));  % Record the differences between their trait
    
    deg_old=nnz(match_PH(idP_old,:)); % number of partners of P_j 
    pr=1/deg_old^eta;  % rewiring probability 
    
    if rand>pr
        idP_z=find(match_PH(:,idH)==0, 1); % find out all plants which the animal unconnected
        if ~isempty(idP_z)
            TD_Hi_z=TD_HP(idH,:); % trait difference between H_i and all its potential plant partners
            TD_new = min(TD_Hi_z(TD_Hi_z>0)); % find out the smallest trait difference
            idP_new = find(TD_Hi_z==TD_new,1); % The ID number of potential plants(P_k) with the best trait match to H_i

            if TD_new < TD_old  % Compare trait matching
                match_PH(idP_old,idH)=0; % disconnect the original interaction
                match_PH(idP_new,idH)=trait_match(Z_P(idP_new),Z_H(idH),sigma); % rewiring, update their trait matching
            end
        end
    end
    
else  % rewiring M->P: pollinator to  plant
    % The following steps are the same as above
    idM=randi(Sm);
    while isempty(find(match_PM(:,idM),1))  
        idM=randi(Sm);
    end
    
    idP_old=index_P2A(idM+Sh);
    TD_old=abs(Z_M(idM)-Z_P(idP_old));  
    
    deg_old=nnz(match_PM(idP_old,:));
    pr=1/deg_old^eta; 
    if rand>pr
        idP_z=find(match_PM(:,idM)==0, 1);
        if ~isempty(idP_z)
            TD_Mi_z=TD_MP(idM,:);
            TD_new = min(TD_Mi_z(TD_Mi_z>0));
            idP_new = find(TD_Mi_z==TD_new,1);
            if TD_new < TD_old
                match_PM(idP_old,idM)=0;
                match_PM(idP_new,idM)=trait_match(Z_P(idP_new),Z_M(idM),sigma);
            end
        end
    end
    
end

end