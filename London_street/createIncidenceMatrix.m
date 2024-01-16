function [E, edgelist] = createIncidenceMatrix(A)
% For a given undirected graph with adj matrix A
% create node-edge incidence matrix E
% Note that self loops are not taken into account here

% compute number of links
links= nnz(A)/2; % number of edges
nodes = length(A);

% find all indices of nodes
[i,j, v] = find(triu(A));
% allocate incidence matrix
% insert +1/-1 in both node rows and column given by edge id
E = sparse(i,1:links,-1,nodes,links);
E = E + sparse(j,1:links,1,nodes,links);

% edgelist in format edge_id, from, to, strength
edgelist =[(1:links)' i j  v];

end
