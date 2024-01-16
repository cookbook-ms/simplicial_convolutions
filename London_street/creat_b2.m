function [b2, triangle_list] = creat_b2(A,edgelist)
% For a given undirected graph with adj matrix A
% create the edge-triangle incidence matrix T
triangle_list = [];
for u = 1:length(A)
    nbs = [edgelist(edgelist(:,2)==u,3);edgelist(edgelist(:,3)==u,2)]; % find the neighbours of each node
    for j = 1:length(nbs)
        for k = j+1:length(nbs)
            v = nbs(j); w = nbs(k); % for each two neighbours, check if they are connected
            nbs_v = [edgelist(edgelist(:,2)==v,3);edgelist(edgelist(:,3)==v,2)]; % find v's neighbours
            nbs_w = [edgelist(edgelist(:,2)==w,3);edgelist(edgelist(:,3)==w,2)]; % find w's neighbours
            if ismember(v,nbs_w) && ismember(w,nbs_v)
                triangle_list = [triangle_list;[u,v,w,1]];
            end
        end
    end
end

% there are repetition in the list
to_be_removed = [];
for i = 1:length(triangle_list)
    if ismember(i,to_be_removed) == 0
        for j = i+1:length(triangle_list)
            if isempty(setdiff(triangle_list(i,1:3), triangle_list(j,1:3)))
                to_be_removed = [to_be_removed;j];
            end
        end
    end
end
triangle_list(to_be_removed,:) = [];
% triangle_list: triangle index, node u, node v, node w, weight = 1
triangle_list = [(1:length(triangle_list))' triangle_list];

rows=[];
cols=[];
vals = [];
% build the incidence matrix B2 (T)
for i = 1:length(triangle_list)
    for j = 1:length(edgelist)
        if triangle_list(i,2:3) == edgelist(j,2:3)
            %T(i,j) = 1;
            rows = [rows;i]; cols=[cols;j]; vals=[vals;1];
        elseif triangle_list(i,3:4) == edgelist(j,2:3)
            %T(i,j) = 1;
            rows = [rows;i]; cols=[cols;j]; vals=[vals;1];
        elseif [triangle_list(i,2),triangle_list(i,4)] == edgelist(j,2:3)
            %T(i,j) = -1;
            rows = [rows;i]; cols=[cols;j]; vals=[vals;-1];
        end
    end
end
T = sparse(rows,cols,vals,length(triangle_list),length(edgelist));
b2 = T';
end
