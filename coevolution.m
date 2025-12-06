% This function is used for the coevolution of traits,
% and records the ID number of the plant that has the least evolutionary effect on each animal.

function [trait_new,idx_P2A]=coevolution(ymea,Sp,Sh,Sm,trait,OptE,omega_ant,omega_mut,omega_env,match_PM,match_PH,epsilon,phi)

trait_new=zeros(Sp+Sh+Sm,1); % new traits for plants, herbivores, and pollinators
idx_P2A=zeros(Sh+Sm,1);      % for each animal, find the plant that has the least evolutionary effect on it

%%
% current traits
Z_P=trait(1:Sp);             % plant
Z_H=trait(Sp+1:Sp+Sh);       % herbivore
Z_M=trait(Sp+Sh+1:Sp+Sh+Sm); % pollinator

% phenotype favoured by the environmental selection
OptE_P=OptE(1:Sp);
OptE_H=OptE(Sp+1:Sp+Sh);
OptE_M=OptE(Sp+Sh+1:Sp+Sh+Sm); 

% evolutionary speed
phi_P=phi(1:Sp);
phi_H=phi(Sp+1:Sp+Sh); 
phi_M=phi(Sp+Sh+1:Sp+Sh+Sm);

% current abundance 
PP=ymea(1:Sp); 
HH=ymea(Sp+1:Sp+Sh);  
MM=ymea(Sp+Sh+1:Sp+Sh+Sm); 

omega_s=1-omega_env;

%%
% calculate phenotype selection 
I_PM=zeros(Sp,Sm);
 for i=1:Sp
    for j=1:Sm
        if match_PM(i,j)>0
            I_PM(i,j)=Z_M(j)-Z_P(i);
        end
    end
 end 
%
I_MP=-I_PM';
%
I_HP=zeros(Sh,Sp);
for i=1:Sh
    for j=1:Sp
        if match_PH(j,i)>0
            I_HP(i,j)=Z_P(j)-Z_H(i);
        end
    end
end
%
I_PH=zeros(Sp,Sh);
for i=1:Sp
    for j=1:Sh
        if (match_PH(i,j)>0 && abs(Z_H(j)-Z_P(i))<epsilon)
            if Z_H(j)>Z_P(i)
                I_PH(i,j)=Z_H(j)-Z_P(i)-epsilon;
            else 
                I_PH(i,j)=Z_H(j)-Z_P(i)+epsilon;
            end
        end
    end
end

%%
% calculate evolutionary effect, 
% and records the ID number of the plant that has the least evolutionary effect on each animal
q_HP=repmat(PP',[Sh,1]).*(match_PH')./(repmat(match_PH'*PP,[1,Sp]));
q_MP=repmat(PP',[Sm,1]).*(match_PM')./(repmat(match_PM'*PP,[1,Sp]));
q_PH=repmat(HH',[Sp,1]).*match_PH./(repmat(match_PH*HH,[1,Sh]));
q_PM=repmat(MM',[Sp,1]).*match_PM./(repmat(match_PM*MM,[1,Sm]));

q_HP(q_HP==0)=NaN;
q_MP(q_MP==0)=NaN;
q_PH(q_PH==0)=NaN;
q_PM(q_PM==0)=NaN;
[min_qHP,idx_qHP]=min(q_HP,[],2);
[min_qMP,idx_qMP]=min(q_MP,[],2);
idx_P2A(1:Sh)=idx_qHP;
idx_P2A(Sh+1:Sh+Sm)=idx_qMP;

q_HP(isnan(q_HP))=0;
q_MP(isnan(q_MP))=0;
q_PH(isnan(q_PH))=0;
q_PM(isnan(q_PM))=0;

%%
% calculate new traits 
Z_P_new=Z_P+phi_P.*(omega_mut*sum(q_PM.*I_PM,2)+omega_ant*sum(q_PH.*I_PH,2)+omega_env*(OptE_P-Z_P));
Z_H_new=Z_H+phi_H.*(omega_s*sum(q_HP.*I_HP,2)+omega_env*(OptE_H-Z_H));
Z_M_new=Z_M+phi_M.*(omega_s*sum(q_MP.*I_MP,2)+omega_env*(OptE_M-Z_M));

trait_new(1:Sp)=Z_P_new;
trait_new(Sp+1:Sp+Sh)=Z_H_new;
trait_new(Sp+Sh+1:Sp+Sh+Sm)=Z_M_new;