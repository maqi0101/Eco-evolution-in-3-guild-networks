% This function is used to generate a random interaction matrix,
% and ensures that no species is isolated.

function matrix_PA=get_matrix_rand(N_row,N_cln,C0)

z_idx=0;
while z_idx==0
    matrix_PA=zeros(N_row,N_cln);

    % Randomly assign interactions based on connection probability C0
    for i=1:N_row
        for j=1:N_cln
            if rand < C0
                matrix_PA(i,j)=1;
            end
        end
    end

    z_cln=all(sum(matrix_PA)~=0);   % check whether any column has an isolated species
    z_row=all(sum(matrix_PA,2)~=0); % Check whether any row has an isolated species

    % If both rows and columns contain no isolated species (z_idx = 1),exit the loop; 
    % otherwise, continue generating the matrix.
    z_idx=z_cln*z_row;  
    
end