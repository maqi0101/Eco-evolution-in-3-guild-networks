% This function is used to describe ecological processes
% systems of differential equations

function dydt = LG_PHM(t,y,Sp,Sh,Sm,r,delta,alpha_PH,gamma_PM,e,h)
%%
dydt=zeros(Sp+Sh+Sm,1); % Sp,Sh,Sm: species number of plant, herbivore, and pollinator
Pi=y(1:Sp);             % plant abundances
Hi=y(Sp+1:Sp+Sh);       % herbivore abundances
Mi=y(Sp+Sh+1:Sp+Sh+Sm); % pollinator abundances

%%
% Interaction relationship matrix
matrix_PH=logical(alpha_PH);
matrix_PM=logical(gamma_PM);

%%
% The antagonistic effect of herbivores on plants
F_ph=matrix_PH.*(Pi*Hi')./(1+h*repmat((matrix_PH'*Pi)',[Sp,1]));

% The mutualistic benefits obtained by pollinators from plants
E1_pm=matrix_PM.*(Pi*Mi')./(1+h*repmat((matrix_PM'*Pi)',[Sp,1]));

% The mutualistic benefits obtained by plants from pollinator
E2_mp=matrix_PM'.*(Mi*Pi')./(1+h*repmat((matrix_PM*Mi)',[Sm,1]));

%%
% Intrinsic growth rate of plants, herbivore, and pollinator
r_P=r(1:Sp);             % plant
r_H=r(Sp+1:Sp+Sh);       % herbivore
r_M=r(Sp+Sh+1:Sp+Sh+Sm); % pollinator

% Self-limiting of plants, herbivore, and pollinator
delt_P=delta(1:Sp);
delt_H=delta(Sp+1:Sp+Sh);
delt_M=delta(Sp+Sh+1:Sp+Sh+Sm);

%%
% Expressions of differential equations
dPi=Pi.*(r_P-delt_P.*Pi)+(sum(gamma_PM'.*E2_mp))'-sum(alpha_PH.*F_ph,2);
dHi=Hi.*(r_H-delt_H.*Hi)+(sum(e*alpha_PH.*F_ph))';
dMi=Mi.*(r_M-delt_M.*Mi)+(sum(gamma_PM.*E1_pm))';

dydt(1:Sp)=dPi;
dydt(Sp+1:Sp+Sh)=dHi;
dydt(Sp+Sh+1:Sp+Sh+Sm)=dMi;

%%
% Temporarily protect species biomass from declining below a near-zero threshold.
% The trigger probability is extremely small and serves as numerical protection.
for i=1:Sp+Sh+Sm
    if y(i)<1e-15 && dydt(i)<0
        dydt(i)=0;
    end
end