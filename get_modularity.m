% This function is used to calculate network Modularity 

function qb = get_modularity(matrix_PA)

% Matlab library "BiMat" is used for analyzing the network structures, 
% which is publicly available via the link: https://bimat.github.io/

fp = Bipartite(matrix_PA);

% Modularity
fp.community = LeadingEigenvector(fp.matrix);  % Newman 2006
fp.community.DoKernighanLinTunning = true;  
fp.community.Detect();
qb = fp.community.Qb;
