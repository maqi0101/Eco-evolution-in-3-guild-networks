% This code is the core program for adaptive eco-evolutionary model, 
% which is used to generate data for ternary diagrams

clc;
clear; 

Sp=30; % species number of plant
Sh=30; % species number of herbivore
Sm=30; % species number of pollinator

r=zeros(Sp+Sh+Sm,1);    % Intrinsic growth rates
delt=zeros(Sp+Sh+Sm,1); % Self-limiting of species
h=0.1;   % Half-saturation constant for antagonistic and mutualistic interactions
e=0.8;   % conversation coefficient of predation 
mu=0.2;  % interaction factor for mutualistic and antagonistic interactions
epsilon=0.3;   % trait barrier 
sigma=0.1;     % controls the sensitivityof trait matching
eta=1;         % parameter for probability of link disconnection

%%
int_ee=0.1:0.1:1;  % range of environmental selection
int_mm=0:0.1:0.9;  % range of mutualistic selection
int_aa=0:0.1:0.9;  % range of antagonistic selection

%%
cM = 4/(Sp+Sm)^0.8;   % mutualistic connectance  
cA = 4/(Sp+Sh)^0.8;   % antagonistic connectance

%% Simulation setup

len=length(int_ee); 
rep=60; % Number of repetitions

% Preallocate arrays to store network and trait properties
DCT=zeros(len,len,rep);   % directionality of traits evolution
FDQ=zeros(len,len,rep);   % functional diversity
SFT=zeros(len,len,rep);   % total functional turnover

Z_cv=zeros(len,len,rep);  % CV of traits
Z_PST=zeros(len,len,rep); % temporal stability of traits
lam=zeros(len,len,rep);   % resilience of networks

ant_N=zeros(len,len,rep); % antagonistic nestedness 
mut_N=zeros(len,len,rep); % mutualistic nestedness
ant_Q=zeros(len,len,rep); % antagonistic modularity
mut_Q=zeros(len,len,rep); % mutualistic modularity

corr_LD=zeros(len,len,rep);  % correlation between L_mut and L_ant
corr_DC=zeros(len,len,rep);  % correlation between D_mut and D_ant
corr_BIO=zeros(len,len,rep); % correlation between B_mut and B_ant
corr_D_B=zeros(len,len,rep); % correlation between (D_mut-D_ant) and (B_mut-B_ant)

%%
counter_max = 100000; % total number of coevolution
m_tep=5000;           % last 5000 iterations for analysis

%% Main simulation loop
for kk=1:rep

    % Random initialization of phenotypes and evolutionary speeds
    OptE=rand(Sp+Sh+Sm,1);            % phenotype favoured by the environmental selection
    trait_org=rand(Sp+Sh+Sm,1);       % original traits of each species
    phi=normrnd(0.5,0.01,Sp+Sh+Sm,1); % evolutionary speed

    % intrinsic growth rate 
    r(1:Sp)=normrnd(1,0.01,Sp,1);             % plants
    r(Sp+1:Sp+Sh)=normrnd(1,0.01,Sh,1);       % herbivore
    r(Sp+Sh+1:Sp+Sh+Sm)=normrnd(1,0.01,Sm,1); % pollinator

    % self-limiting of plants, herbivore, and pollinator
    delt(1:Sp)=normrnd(1,0.01,Sp,1);
    delt(Sp+1:Sp+Sh)=normrnd(1,0.01,Sh,1);
    delt(Sp+Sh+1:Sp+Sh+Sm)=normrnd(1,0.01,Sm,1);
    
    matrix_PH_origa=get_matrix_rand(Sp,Sh,cA); % original antagonistic subnetwork（0-1 matrix）
    matrix_PM_origa=get_matrix_rand(Sp,Sm,cM); % original mutualistic subnetwork（0-1 matrix）

    for ii=1:len
        for jj=1:ii % iterate over ii and jj: 55 specific combinations of selection pressures        
            omega_ant=int_aa(ii-jj+1);  % Weight of antagonistic selection
            omega_mut=int_mm(jj);       % Weight of mutualistic selection
            omega_env=int_ee(len-ii+1); % Weight of environmental selection
            
            trait=trait_org; % original traits

            % trait matching between plant and herbivore
            match_PH=zeros(Sp,Sh); 
            for i=1:Sp
                for j=1:Sh
                    match_PH(i,j)=matrix_PH_origa(i,j)*trait_match(trait(i),trait(Sp+j),sigma);
                end
            end

            % trait matching between plant and pollinator
            match_PM=zeros(Sp,Sm); 
            for i=1:Sp
                for j=1:Sm
                    match_PM(i,j)=matrix_PM_origa(i,j)*trait_match(trait(i),trait(Sp+Sh+j),sigma);
                end
            end
            
            gamma_PM=mu*match_PM; % strength of mutualistic interactions 
            alpha_PH=mu*match_PH; % strength of antagonistic interactions
            
            mu_N0 = 0.1;
            N0 = mu_N0*ones(Sp+Sh+Sm,1);  % initial abundance
            to = 0;    % initial time
            tf = 50;   % integration time within each iteration 
            counter=0; % current number of iterations
            
            fdq_new=0; % record current functional diversity
            fdq_old=0; % functional diversity in the previous iteration

            % several properties for ecological networks of the last 5000 iterations
            FDQ_ijk=zeros(1,m_tep);  % functional diversity
            SFT_ijk=zeros(1,m_tep);  % functional turnover
            lam_ijk=zeros(1,m_tep);  % resilience of networks

            ant_Q_ijk=zeros(1,m_tep); % antagonistic modularity
            ant_N_ijk=zeros(1,m_tep); % antagonistic weight nestedness
            mut_Q_ijk=zeros(1,m_tep); % mutualistic modularity
            mut_N_ijk=zeros(1,m_tep); % mutualistic weight nestedness
            
            corr_LD_ijk=zeros(1,m_tep);  % correlation between L_mut and L_mut
            corr_DC_ijk=zeros(1,m_tep);  % correlation between D_mut and D_mut
            corr_BIO_ijk=zeros(1,m_tep); % correlation between B_mut and B_mut
            corr_D_B_ijk=zeros(1,m_tep); % correlation between (D_mut-D_mut) and (B_mut-B_mut)


            ZZ_ijk=zeros(Sp+Sh+Sm,m_tep); % species traits during the last 5000 iterations
            ZZ_chg=zeros(Sp+Sh+Sm,counter_max); % trait change across all iterations
            
            
            while counter<=counter_max
                
                if counter>0
                    % coevolution; and for each animal, find the plant that has the least evolutionary effect on it
                    [trait_new,idx_P2A]=coevolution(N0,Sp,Sh,Sm,trait,OptE,omega_ant,omega_mut,omega_env,match_PM,match_PH,epsilon,phi);
                    
                    ZZ_chg(:,counter)=trait_new-trait; % traits change
                    trait=trait_new; % update traits

                    % update trait matching, and rewiring one link
                    [match_PH,match_PM]=updata_interaction(Sp,Sh,Sm,trait,idx_P2A,match_PM,match_PH,sigma,eta);
                    
                    gamma_PM=mu*match_PM; % update the strength of interspecific interactions (mutualism)
                    alpha_PH=mu*match_PH; % update the strength of interspecific interactions (antagonism)
                end
                
                % integration (ecological processes)
                [t y] = ode45(@(t,y)LG_PHM(t,y,Sp,Sh,Sm,r,delt,alpha_PH,gamma_PM,e,h),[to to+tf],N0);
                to = t(end);
                N0=y(end,:)';  % set initial condition for next interval
                
                %
                if counter>=counter_max-m_tep
                    fdq_new= get_FDQ(trait,N0); % current functional diversity
                end
                
                %
                if counter>counter_max-m_tep

                    % functional diversity and functional turnover
                    FDQ_ijk(counter+m_tep-counter_max)=fdq_new;
                    SFT_ijk(counter+m_tep-counter_max)=abs(fdq_new-fdq_old);
                    
                    % resilience
                    jacob_mat=get_jacmat(N0,Sp,Sh,Sm,r,delt,alpha_PH,gamma_PM,e,h);
                    lam_ijk(counter+m_tep-counter_max)=-max(real(eig(jacob_mat)));
                    
                    % calculate modularity and nestedness
                    ant_Q_ijk(counter+m_tep-counter_max)=get_modularity(alpha_PH);
                    mut_Q_ijk(counter+m_tep-counter_max)=get_modularity(gamma_PM);
                    ant_N_ijk(counter+m_tep-counter_max)=get_wNoda(alpha_PH);
                    mut_N_ijk(counter+m_tep-counter_max)=get_wNoda(gamma_PM);
                    
                    ZZ_ijk(:,counter+m_tep-counter_max)=trait; % species traits
                    P_abud=N0(1:Sp);  % abundance of plants 
                    H_abud=N0(Sp+1:Sp+Sh);   % abundance of herbivores
                    M_abud=N0(Sp+Sh+1:Sp+Sh+Sm); % abundance of pollinator
                    
                    % interaction relationship; plant’ degree centrality
                    mat_pm=logical(gamma_PM);
                    mat_ph=logical(alpha_PH);
                    dc_mut=sum(mat_pm,2)'/Sm;
                    dc_ant=sum(mat_ph,2)'/Sh;
                    
                    % energy flow situation
                    mut_mat=(diag(M_abud)*gamma_PM')./(1+h*repmat((mat_pm*M_abud)',[Sm,1]));
                    ant_mat=(alpha_PH*diag(H_abud))./(1+h*repmat((mat_ph'*P_abud)',[Sp,1]));
                    bio_mut=sum(mut_mat);
                    bio_ant=sum(ant_mat,2)';
                    
                    % link diversity
                    gamma_PM_w=gamma_PM./repmat(sum(gamma_PM,2),[1,Sm]);
                    alpha_PH_w=alpha_PH./repmat(sum(alpha_PH,2),[1,Sh]);
                    dh_mut=gamma_PM_w.*log(gamma_PM_w);
                    dh_ant=alpha_PH_w.*log(alpha_PH_w);
                    dh_mut(isnan(dh_mut))=0;
                    dh_ant(isnan(dh_ant))=0;
                    LD_mut=-sum(dh_mut,2)';
                    LD_ant=-sum(dh_ant,2)';
                    
                    % Calculate correlation
                    corr_LD_ijk=corr(LD_mut',LD_ant','Type','Spearman');
                    corr_DC_ijk=corr(dc_mut',dc_ant','Type','Spearman');
                    corr_BIO_ijk=corr(bio_mut',bio_ant','Type','Spearman');
                    corr_D_B_ijk=corr((dc_mut-dc_ant)',(bio_mut-bio_ant)','Type','Spearman');
  
                end
                
                fdq_old=fdq_new;
                fprintf('\nrepeat=%d, row = %d, list = %d,  count=%d\n\n', kk, ii, jj, counter);
                counter = counter+1;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % the mean of the last 5000 iterations
            lam(ii,jj,kk)=mean(lam_ijk);
            FDQ(ii,jj,kk)=mean(FDQ_ijk);
            ant_N(ii,jj,kk)=mean(ant_N_ijk);
            mut_N(ii,jj,kk)=mean(mut_N_ijk);
            ant_Q(ii,jj,kk)=mean(ant_Q_ijk);
            mut_Q(ii,jj,kk)=mean(mut_Q_ijk);
            corr_LD(ii,jj,kk)=mean(corr_LD_ijk);
            corr_DC(ii,jj,kk)=mean(corr_DC_ijk);
            corr_BIO(ii,jj,kk)=mean(corr_BIO_ijk);
            corr_D_B(ii,jj,kk)=mean(corr_D_B_ijk);

            % total functional turnover
            SFT(ii,jj,kk)=sum(SFT_ijk);

            % directionality of traits evolution
            DCT(ii,jj,kk)=sum(abs(trait-trait_org)./sum(abs(ZZ_chg),2))/(Sp+Sh+Sm);
            DCT(isnan(DCT))=0;
            
            % CV of the total amount of normalized traits 
            z_max=max(max(ZZ_ijk));
            z_min=min(min(ZZ_ijk));
            z_nor=(ZZ_ijk-z_min)/(z_max-z_min);
            Z_cv(ii,jj,kk)=std(sum(z_nor))/mean(sum(z_nor));
            
            % temporal stability of traits
            if Z_cv(ii,jj,kk)<0.005
                Z_PST(ii,jj,kk)=1;
            end
             
        end                                                                                                                                                                                                                                                       
    end
end

%%
% the mean of 60 repetitions
DCT_mean=mean(DCT,3);
FDQ_mean=mean(FDQ,3);
SFT_mean=mean(SFT,3);

lam_mean=mean(lam,3);
Z_cv_mean=mean(Z_cv,3);
Z_PST_mean=mean(Z_PST,3);

ant_N_mean=mean(ant_N,3);
mut_N_mean=mean(mut_N,3);
ant_Q_mean=mean(ant_Q,3);
mut_Q_mean=mean(mut_Q,3);

corr_LD_mean=mean(corr_LD,3);
corr_DC_mean=mean(corr_DC,3);
corr_BIO_mean=mean(corr_BIO,3);
corr_D_B_mean=mean(corr_D_B,3);

%%
% Arrange 55 combinations in a row according to rules and export them as txt
Arr_DCT=zeros(1,55);
Arr_FDQ=zeros(1,55);
Arr_SFT=zeros(1,55);
Arr_Res=zeros(1,55);
Arr_Zcv=zeros(1,55);
Arr_Zpst=zeros(1,55);
Arr_Nant=zeros(1,55);
Arr_Nmut=zeros(1,55);
Arr_Qant=zeros(1,55);
Arr_Qmut=zeros(1,55);
Arr_Cor_ld=zeros(1,55);
Arr_Cor_dc=zeros(1,55);
Arr_Cor_bio=zeros(1,55);
Arr_Cor_d_b=zeros(1,55);
%
sot=0;
for jj=1:10
    for ii=10:-1:jj
        sot=sot+1;
        Arr_DCT(sot)=DCT_mean(ii,jj);
        Arr_FDQ(sot)=FDQ_mean(ii,jj);
        Arr_SFT(sot)=SFT_mean(ii,jj);
        Arr_Res(sot)=lam_mean(ii,jj);
        Arr_Zcv(sot)=Z_cv_mean(ii,jj);
        Arr_Zpst(sot)=Z_PST_mean(ii,jj);
        Arr_Nant(sot)=ant_N_mean(ii,jj);
        Arr_Nmut(sot)=mut_N_mean(ii,jj);
        Arr_Qant(sot)=ant_Q_mean(ii,jj);
        Arr_Qmut(sot)=mut_Q_mean(ii,jj);
        Arr_Cor_ld(sot)=corr_LD_mean(ii,jj);
        Arr_Cor_dc(sot)=corr_DC_mean(ii,jj);
        Arr_Cor_bio(sot)=corr_BIO_mean(ii,jj);
        Arr_Cor_d_b(sot)=corr_D_B_mean(ii,jj);

    end
end

%%
% Save the corresponding indicators as a txt file, 
% and then draw it as a ternary diagrams using a python program called Draw_Ternary.py

writematrix(Arr_DCT, 'Directionality.txt', 'Delimiter', ',');
writematrix(Arr_FDQ, 'FDQ.txt', 'Delimiter', ',');
writematrix(Arr_SFT, 'SFT.txt', 'Delimiter', ',');
writematrix(Arr_Res, 'Resilience.txt', 'Delimiter', ',');
writematrix(Arr_Zcv, 'Trait_cv.txt', 'Delimiter', ',');
writematrix(Arr_Zpst, 'Trait_pst.txt', 'Delimiter', ',');
writematrix(Arr_Nant, 'Nant.txt', 'Delimiter', ',');
writematrix(Arr_Nmut, 'Nmut.txt', 'Delimiter', ',');
writematrix(Arr_Qant, 'Qant.txt', 'Delimiter', ',');
writematrix(Arr_Qmut, 'Qmut.txt', 'Delimiter', ',');
writematrix(Arr_Cor_ld, 'Corr_LD.txt', 'Delimiter', ',');
writematrix(Arr_Cor_dc, 'Corr_DC.txt', 'Delimiter', ',');
writematrix(Arr_Cor_bio, 'Corr_BIO.txt', 'Delimiter', ',');
writematrix(Arr_Cor_d_b, 'Corr_D_B.txt', 'Delimiter', ',');
 
