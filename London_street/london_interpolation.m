
%% NOT WORKING
%
% london network interpolation

clc;
clear;
close all;
rng(1223)
addpath('/Users/maosheng/Documents/chebfun')
data = importdata('LondonEdges.csv');
coordinate = importdata('LondonNodes.csv');
x_coord = coordinate.data(:,2);
y_coord = coordinate.data(:,3);
from = str2double(data.textdata(2:end,1))';
to = str2double(data.textdata(2:end,2))';
num_nodes = max([from, to]);

% num_edges = length(from);
% B1 = sparse(from,1:num_edges,-1,num_nodes,num_edges) ...
%     + sparse(to,1:num_edges,1,num_nodes,num_edges);
% some edges are encoded in both directions -- we need to eliminate them
% for this example, so we construct the adjacency matrix first, and from
% there get the incidence matrix.
A = sparse(from,to,1,num_nodes,num_nodes);
A = A+A';
A(A>1) = 1;

[B1, edgelist] = createIncidenceMatrix(A);
from = edgelist(:,2)';
to = edgelist(:,3)';
num_edges= size(B1,2);
B1 = full(B1);
% plot the London street network
p1 = figure;
plt_network(B1,ones(num_edges,1),x_coord,y_coord);
axis tight

% creat the edge-triangle incidence matrix
[B2, triangle_list] = creat_b2(A,edgelist);
B2 = full(B2);
num_triangles = size(B2,2);

% Hodge Laplacian
L1l = B1'*B1; L1u = B2*B2';
L1 = L1l + L1u;

% projectors;
P_g = B1'*pinv(B1*B1')*B1;
P_c = B2*pinv(B2'*B2)*B2';
P_h = eye(num_edges) - L1*pinv(L1);
% edge Laplacian
Le = L1l;
% Line graph (undirected) adjacency)
AL = abs(Le-2*eye(size(Le)));
LL = diag(sum(AL,2)) - AL;
% graph Laplacian
L0 = B1*B1';
% triangle Laplacian
L2 = B2'*B2;

% spectrum
[u1l, lam_l]  = eig(L1l); eig_L1l = diag(lam_l); u_g = u1l(:,50:end);
[u1, lam]= eig(L1); eig_L1 = diag(lam); u_h = u1(:,1:37);
[u1u, lam_u]= eig(L1u); eig_L1u = diag(lam_u); u_c = u1u(:,119:end);
lam_g = eig_L1l(50:end); lam_c = eig_L1u(119:end);
% note that u_c don't necessarily have the same ordering as in u1;
%% generating mask
ratio = 0.75;
M = floor(num_edges*ratio); % the number of nonzero nodes
mask = zeros(num_edges,1);
mask(randperm(numel(mask), M)) = 1;

%% data drive method to learn the interpolation filter
num_training = 100;
num_test = 50;
out = randn(num_edges,num_training);
in = mask.*out;

out = reshape(out, [], 1);

% data driven to learn the interpolation filter
l1 = 1:10; l2 = 1:10;
sys_mat1 = [];
for i = 1:num_training
    aa = [in(:,i)];
    for l1i = l1
        aa = [aa L1l^l1i*in(:,i)];
    end
    sys_mat1 = [sys_mat1;aa];
end
sys_mat2 = [];
for i = 1:num_training
    aa = [];
    for l2i = l2
        aa = [aa L1u^l2i*in(:,i)];
    end
    sys_mat2 = [sys_mat2;aa];
end
sys_mat = [sys_mat1 sys_mat2];

% obtain the filter coefficient
h = pinv(sys_mat)*out;
h = inv(sys_mat'*sys_mat + 0.1*eye(size(sys_mat'*sys_mat)))*sys_mat'*out;
H = h(1)*eye(num_edges);
for l = l1
    H = H + h(1+l)*L1l^l;
end
for l = l2
    H = H + h(1+l1(end)+l)*L1u^l;
end

 
err = 0;
for i = 1:num_test
    test_out = randn(num_edges,1);
    test_in = mask.*test_out;
    test_out_approx = H*test_in;
    err_original(i) = norm(test_out-test_in)/norm(test_out); 
    err(i) = norm(test_out - test_out_approx)/norm(test_out);
end
mean(err)
mean(err_original)