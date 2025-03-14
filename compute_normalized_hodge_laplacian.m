function [L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2)
%This is a small function to compute the normalized hodge laplacian based
%on the 1- and 2- incidence matrices of a simplicial complex
% the method is based on the paper: rankdom walks and the normalized hodge
% 1-Laplacian by Michael T. Schaub
[num_nodes, num_edges] = size(B1);
num_tri = size(B2,2);
D2 = max(diag(abs(B2)*ones(num_tri,1)),eye(num_edges));
D1 = 2*diag(abs(B1)*D2*ones(num_edges,1));
D3 = 1/3*eye(num_tri);
L1_n = D2*B1'/D1*B1 + B2*D3*B2'/D2;
end

