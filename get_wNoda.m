% This function is used to calculate weighted Nestedness (WNODA; Pinheiro 2019)

function wNoda = get_wNoda(mat_PA) 

% mat_PA is the interaction strength matrix

[SSp,SSa]=size(mat_PA); % number of rows and columns

% Sort rows and columns by decreasing marginal totals
num_row=sum(mat_PA,2)'; % total marginal number of rows
num_cln=sum(mat_PA,1);  % total marginal number of columns

[val1,idx1]=sort(num_row,'descend'); % sort row totals in descending order
[val2,idx2]=sort(num_cln,'descend'); % sort column totals in descending order

mat_PA1=mat_PA(idx1,:);     % matrix sorted by row marginal totals
Ord_mat_PA=mat_PA1(:,idx2); % sorted by decreasing marginal totals of columns

% nestedness between rows
wNoda_row=0; 
for ii=1:SSp-1
    for jj=ii+1:SSp
        row_i=Ord_mat_PA(ii,:); % row i
        row_j=Ord_mat_PA(jj,:); % row j
        
        nub_i=val1(ii);      % marginal total of row i
        nub_j=val1(jj);      % marginal totals of row j
        deg_j=sum(row_j~=0); % number of non-zero elements in row j  
        
        % calculate weighted pairwise overlap
        if nub_i>nub_j && nub_j>0 
            small_ij=zeros(1,SSa);
            for kk=1:SSa
                if row_j(kk)>0 && row_j(kk)<row_i(kk)
                    small_ij(kk)=1;
                end
            end
            k_ij=100*sum(small_ij)/deg_j;
        else
            k_ij=0;
        end
        wNoda_row=wNoda_row+k_ij;
    end
end

% nestedness between columns
% the following steps are the same as above
wNoda_cln=0; 
for ii=1:SSa-1
    for jj=ii+1:SSa
        cln_i=Ord_mat_PA(:,ii);
        cln_j=Ord_mat_PA(:,jj);

        nub_i=val2(ii); 
        nub_j=val2(jj); 
        deg_j=sum(cln_j~=0); 
        
        if nub_i>nub_j && nub_j>0
            small_ij=zeros(1,SSp);
            for kk=1:SSp
                if cln_j(kk)>0 && cln_j(kk)<cln_i(kk)
                    small_ij(kk)=1;
                end
            end
            k_ij=100*sum(small_ij)/deg_j;
        else
            k_ij=0;
        end
        wNoda_cln=wNoda_cln+k_ij;
    end
end

% Overall matrix nestedness (row-based + column-based components)
wNoda=(wNoda_row+wNoda_cln)/(SSp*(SSp-1)/2+SSa*(SSa-1)/2);