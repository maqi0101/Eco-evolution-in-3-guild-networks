% This function is used to get the jacobian matrix

function jac_mat = get_jacmat(ymea,Sp,Sh,Sm,r,delt,alpha_PH,gamma_PM,e,h)
%%
% abundance of each species at ecological equilibrium
PP=ymea(1:Sp);             % plant
HH=ymea(Sp+1:Sp+Sh);       % herbivore
MM=ymea(Sp+Sh+1:Sp+Sh+Sm); % pollinator

% intrinsic growth rates of plants, herbivore, and pollinator
r_P=r(1:Sp);
r_H=r(Sp+1:Sp+Sh);
r_M=r(Sp+Sh+1:Sp+Sh+Sm);

% self-limiting of plants, herbivore, and pollinator
delt_P=delt(1:Sp);
delt_H=delt(Sp+1:Sp+Sh);
delt_M=delt(Sp+Sh+1:Sp+Sh+Sm);

% interaction relationship matrix, 0-1 matrix 
mat_PH=logical(alpha_PH);
mat_PM=logical(gamma_PM);
mat_MP=mat_PM';
mat_HP=mat_PH';

%%
% dP/dP
dif_PP=zeros(Sp,Sp);
for i=1:Sp
    for j=1:Sp
        if (i==j)
            dif_PP(i,j)=PP(i)*(-delt_P(i)+h*alpha_PH(i,:).*mat_PH(j,:)./((1+h*mat_HP*PP)'.^2)*HH);
        else
            dif_PP(i,j)=PP(i)*(h*alpha_PH(i,:).*mat_PH(j,:)./((1+h*mat_HP*PP)'.^2)*HH);
        end
    end
end
 
% dP/dH
dif_PH=zeros(Sp,Sh);
for i=1:Sp
    for j=1:Sh
        dif_PH(i,j)=-PP(i)*alpha_PH(i,j)/(1+h*mat_HP(j,:)*PP);
    end
end
 
% dP/dM
dif_PM=zeros(Sp,Sm);
for i=1:Sp
    for j=1:Sm
        dif_PM(i,j)=PP(i)/(1+h*mat_PM(i,:)*MM)*(gamma_PM(i,j)+mat_PM(i,j)*h*(r_P(i)-delt_P(i)*PP(i)-alpha_PH(i,:)./((1+h*mat_HP*PP)')*HH));
    end
end
 
% dH/dP
dif_HP=zeros(Sh,Sp);
for i=1:Sh
    for j=1:Sp
        dif_HP(i,j)=HH(i)*(e*alpha_PH(j,i)+mat_HP(i,j)*h*(r_H(i)-delt_H(i)*HH(i)))/(1+h*mat_HP(i,:)*PP);
    end
end
 
% dH/dH
dif_HH=zeros(Sh,Sh);
for i=1:Sh
    for j=1:Sh
        if (i==j)
            dif_HH(i,j)=-delt_H(i)*HH(i);
        else
            dif_HH(i,j)=0;
        end
    end
end
 
% dH/dM
dif_HM=zeros(Sh,Sm);
 
% dM/dP
dif_MP=zeros(Sm,Sp);
for i=1:Sm
    for j=1:Sp
        dif_MP(i,j)=MM(i)*(gamma_PM(j,i)+mat_MP(i,j)*h*(r_M(i)-delt_M(i)*MM(i)))/(1+h*mat_MP(i,:)*PP);
    end
end

% dM/dH
dif_MH=zeros(Sm,Sh);

% dM/dM
dif_MM=zeros(Sm,Sm);
for i=1:Sm
    for j=1:Sm
        if (i==j)
            dif_MM(i,j)=-delt_M(i)*MM(i);
        else
            dif_MM(i,j)=0;
        end
    end
end

%%
% jacobian matrix
jac_mat=[dif_PP,dif_PH,dif_PM;
         dif_HP,dif_HH,dif_HM;
         dif_MP,dif_MH,dif_MM];
