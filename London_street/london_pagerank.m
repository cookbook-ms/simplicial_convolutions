% London network pagerank on edges study
clc;
clear;
close all;
rng(1223)
addpath('/Users/maosheng/Documents/chebfun')
addpath('/Users/maosheng/Documents/subtightplot')
addpath('/Users/maosheng/OneDrive - Delft University of Technology/Topological filter/codes from 15 update/experiments')
data = importdata('LondonEdges.csv');
coordinate = importdata('LondonNodes.csv');
x_coord = coordinate.data(:,2);
y_coord = coordinate.data(:,3);
from = str2double(data.textdata(2:end,1))';
to = str2double(data.textdata(2:end,2))';
num_nodes = max([from, to]);

A = sparse(from,to,1,num_nodes,num_nodes);
A = A+A';
A(A>1) = 1;

[B1, edgelist] = createIncidenceMatrix(A);
from = edgelist(:,2)';
to = edgelist(:,3)';
num_edges= size(B1,2);
B1 = full(B1);
 
% creat the edge-triangle incidence matrix
[B2, triangle_list] = creat_b2(A,edgelist);
B2 = full(B2);
B2t = B2';
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

%% PageRank to study the importance of some edges
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% build the pagerank operator
k = 0.01;
H_P = inv(k*eye(num_edges) + L1_n);
% eigendecomposition of L1_n
[U_n,Lam_n,V_n] = eig(L1_n); % U_n is the right eigenvector; V_n' is the left one
Lam_n = diag(Lam_n);
Lam_n(Lam_n(:)<1e-3) = 0;
% approximation design
lam_n = uniquetol(sort(Lam_n),0.1);
% universal design
lam_n = linspace(0,1,201)';
% let us use the normalized Laplacian to formulate a filter to model the
% pagerank operator
% frequency response
h_p = 1./(k+lam_n);
% filter order is 10
for L = 1:10 % filter length ind_norm_h
    Phi_n = lam_n.^(0:1:L-1);
    alpha = pinv(Phi_n) * h_p;
    % build the filter
    H_pr = zeros(num_edges);
    for l = 1:L
        H_pr = H_pr + alpha(l)*L1_n^(l-1);
    end
    err(L) = norm(H_pr - H_P,'fro')/norm(H_P,'fro');
end

% continuous version of the frequency response 
h_p_c = @(x) (1./(k+x));

% true frequency response 
H_P_fr = h_p_c(Lam_n);
H_P_fr = diag(U_n'*H_P*U_n);

% perform the chebyshev design

% obtain the chebyshev approx and coefficients
h_p_cheb = chebfun(h_p_c,[0,1]); % note that the maximal eigenvalue is 1 in normalized laplacian
cheb_coeff = chebcoeffs(h_p_cheb);
alpha_g = 1/2;
% the truncated order
K_trnc = 2:4:length(cheb_coeff);
for i = 1:length(K_trnc)
    H_p_cheb_approx_out(i,:,:) = filter_cheb_approx(L1_n,cheb_coeff,alpha_g,K_trnc(i));
end
 
for i = 1:length(K_trnc)
    h_p_cheb_approx_out = diag(U_n'*squeeze(H_p_cheb_approx_out(i,:,:))*U_n);
    % compute the error w.r.t. the true gradient extraction filter
    % err_cheb(i) = norm(h_p_cheb_approx_out - H_P_fr)/norm(H_P_fr);
    % we could consider the normalized frobenius norm or the 2 norm, the
    % spectral norm
    %err_filter_cheb(i) = norm(squeeze(H_p_cheb_approx_out(i,:,:))-H_P,'fro');
    err_filter_cheb(i) = norm(squeeze(H_p_cheb_approx_out(i,:,:))-H_P,2);
    err_cheb_bound(i) = max(abs(H_P_fr-h_p_cheb_approx_out));
end
% find the chebyshev orders which lead to an error smaller than 1e-3;
K_proper = find(err_filter_cheb(:)<1e-3); K_chosen = K_proper(1);

% identity term coefficient to confirm the correct calculation
sum(cheb_coeff(1:2:end))-sum(cheb_coeff(2:2:end));

%% plot to compare
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.08], [0.1 0.03], [0.08 0.01]);
if ~make_it_tight,  clear subplot;  end

figure;
subplot(2,2,1)
width_ratio = 1;
plt_network(B1,ones(num_edges,1),x_coord,y_coord,width_ratio,0);
axis tight

subplot(2,2,2);
plot(lam_n,h_p,'LineWidth',2); hold on;
% universal designed frequency response or approximatedly designed one
h_pr_universal = diag(U_n'*H_pr*U_n);
scatter(Lam_n, h_pr_universal,50,[0.6350 0.0780 0.1840]); hold on;
% let us plot the one by chebyshev approximation with the first K leading
% to an error smaller than 1e-3
scatter(Lam_n, diag(U_n'*squeeze(H_p_cheb_approx_out(K_chosen,:,:))*U_n),...
    50,[0.4660 0.6740 0.1880], 's');
xlim([0 1]);
ylim([0 100]); 
grid on;
set(gca,'fontsize',14)
set(gca, 'YScale', 'log')
legend(['Continuous desired ','$h_{PR}(\lambda)$'],'Universal design with 200 points',...
    'Chebyshev with K=62',...
    'Interpreter','latex')
set(gca,'xticklabel',{[]})
% plot the chebyshev approximation error as the order increases

subplot(2,2,3);
plot(K_trnc, err_filter_cheb,'LineWidth',2); hold on;
%plot(K_trnc, err_cheb_bound,'--','LineWidth',3); hold on;
scatter(K_trnc(K_chosen),err_filter_cheb(K_chosen),100,'rs','filled')
set(gca, 'YScale', 'log')
grid on;
set(gca,'fontsize',14)
legend('$$||\mathbf{H}_{PR}-\mathbf{H}_{1,c}||_2$$','error of order 62','Interpreter','latex');
xlabel('Chebyshev order')
ylim([min(err_filter_cheb) max(err_filter_cheb)])

subplot(2,2,4);
lam_conti = 0:0.001:1; 
in = (lam_conti-alpha_g)./alpha_g;

for i = 6:8:24
    %h_p_cheb_plot = diag(U_n'*squeeze(H_p_cheb_approx_out(i,:,:))*U_n);
   plot(lam_conti,abs(h_p_c(lam_conti)-cheb_approx(in,cheb_coeff(1:K_trnc(i)))),'LineWidth',2)
    %scatter(Lam_n,abs(h_p-h_p_cheb_plot),'LineWidth',2);
    hold on;
end
set(gca, 'YScale', 'log')
grid on;
set(gca,'fontsize',14)
xlabel('Frequency')
ylabel('$|h_{PR}(\lambda) - h_{1,c}(\lambda)|$','Interpreter','latex')
 legend(['order ',num2str(K_trnc(6))], ['order ',num2str(K_trnc(14))], ...
     ['order ',num2str(K_trnc(22))], 'interpreter','latex')
%% pagerank 1 -- study several major edges
% especially the edges with large upper degree
% deg = diag(L1u);
% % pos tells us which edges have large degree，we consider the first three
% % largest
% [deg_max,pos] = maxk(deg,15);
% can we choose the ones with large number of neighbours to study

% below we use the chebyshev approximation to perform pagerank study
H_pr = squeeze(H_p_cheb_approx_out(K_chosen,:,:)); 

% BASICALLY put the edge index wanted to be studied in POS, then gives the
% most influenced edge indices
num_most = 5; % the top 5 most influenced
for i = 1:num_edges
    num_nbr_l(i) = nnz(L1l(:,i)); % number of neighbors + 1
end
[num_nbr_l_max,pos_nbr_l] = maxk(num_nbr_l,num_most);

for i = 1:num_edges
    num_nbr_u(i) = nnz(L1u(:,i));  % number of neighbors + 1
end
[num_nbr_u_max,pos_nbr_u] = maxk(num_nbr_u,num_most);

% choose two edges to study their influences
pos = [30;42];
% pos = 30;

% used for page ranking
for i = 1:length(pos)
    ind_vec = zeros(num_edges,1);
    ind_vec(pos(i)) = 1;
    % the pagerank vector
    pg_vec(:,i) = H_P*ind_vec;
    % the filtering method
    pg_vec_H(:,i) = H_pr*ind_vec;
    err_pr_vec(L,i) = norm(pg_vec - pg_vec_H)/norm(pg_vec);
    % look for the indices with a large pagerank value larger than 5
    [pg_vec_max5(:,i), ind_pg(:,i)] = maxk(pg_vec(:,i),num_most,...
        'ComparisonMethod','abs');
    [pg_vec_H_max5(:,i), ind_pg_H(:,i)] = maxk(pg_vec(:,i),num_most,...
        'ComparisonMethod','abs');
end

% find the connected nodes of these most influenced edges
for i = 1:length(pos)
    for j = 1:num_most
        ind_pg_node1(j,i) = find(B1(:,ind_pg(j,i))==1);
        ind_pg_node2(j,i) = find(B1(:,ind_pg(j,i))==-1);
    end
end

%% make figures

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.08], [0.00 0.03], [0.00 0.00]);
if ~make_it_tight,  clear subplot;  end
% set 5 colors for 5 most influenced edges
figure;
subplot(1,1,1)
% if more wanted to be shown, then set more colors
hex = ['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30'];
color = string(hex);
for i = 1:length(pos)
    label = string(ind_pg(:,i));
    plt_network_page_rank(B1,ind_pg(:,i),color,label,x_coord,y_coord);
    hold on;
end
 
map = sscanf(flip(hex)','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
rgbplot(map)
hold on
colormap(map)
colorbar('Ticks',[])
axis tight

%% pagerank 2 -- study the influences on the subspace of all edges
[num_nodes, num_edges] = size(B1);
num_tri = size(B2,2);
D2 = max(diag(abs(B2)*ones(num_tri,1)),eye(num_edges));
D1 = 2*diag(abs(B1)*D2*ones(num_edges,1));
D3 = 1/3*eye(num_tri);
L1_n = D2*B1'/D1*B1 + B2*D3*B2'/D2;
% compute the normalized hodge laplacian decomposition
L1_n_s = D2^(0.5)*B1'/D1*B1*D2^(0.5) + D2^(-0.5)*B2*D3*B2'*D2^(-0.5);
issymmetric(L1_n_s)

[U_n,Lam_n,V_n] = eig(L1_n_s); Lam_n = diag(Lam_n);
Lam_n(Lam_n(:)<1e-3) = 0;
% the harmonic space
U_H_n = U_n(:,1:37);

[Ul_n,Lam_l_n] = eig(D2^(0.5)*B1'*(D2^(0.5)*B1')'); Lam_l_n = diag(Lam_l_n);
Lam_l_n(Lam_l_n(:)<1e-3) = 0;
% the gradient space
U_G_n = Ul_n(:,50:end);

[Uu_n,Lam_u_n] = eig(D2^(-0.5)*B2*(D2^(-0.5)*B2)'); Lam_u_n = diag(Lam_u_n);
Lam_u_n(Lam_u_n(:)<1e-3) = 0;
% the curl space
U_C_n = Uu_n(:,119:end);
% for each edge study its influence
for i = 1:num_edges
    ind_vec = zeros(num_edges,1);
    ind_vec(i) = 1;
    % the pagerank vector
    pg_vec(:,i) = H_P*ind_vec;
    % the filtering method
    pg_vec_H(:,i) = H_pr*ind_vec;
    % evaluate in an aggregated fashion through 2-norm
    pg_norm(i) = norm(pg_vec(:,i));
    pg_H_norm(i) = norm(pg_vec_H(:,i));
    %err_pr_vec(i) = norm(pg_vec - pg_vec_H)/norm(pg_vec);
    % decompose and extract the gradient component, for its norm
    % we use the apprxoimately designed gradient preserving filter
    pg_vec_g(:,i) = U_G_n'*pg_vec(:,i);
    pg_vec_H_g(:,i) = U_G_n'*pg_vec_H(:,i);
    % evaluate the curl norm
    pg_vec_c(:,i) = U_C_n'*pg_vec(:,i);
    pg_vec_H_c(:,i) = U_C_n'*pg_vec_H(:,i);
    % evaluate the harmonic norm
    pg_vec_h(:,i) = U_H_n'*pg_vec(:,i);
    pg_vec_H_h(:,i) = U_H_n'*pg_vec_H(:,i);

    % the absolute
    pg_norm_g(i) = norm(pg_vec_g(:,i));
    pg_H_norm_g(i) = norm(pg_vec_H_g(:,i));
    % the absolute
    pg_norm_c(i) = norm(pg_vec_c(:,i));
    pg_H_norm_c(i) = norm(pg_vec_H_c(:,i));
    % the absolute
    pg_norm_h(i) = norm(pg_vec_h(:,i));
    pg_H_norm_h(i) = norm(pg_vec_H_h(:,i));
end

% we should evaluate the normalized norm, i.e., wrt the total norm
% evaluate the gradient norm

for i = 1:num_edges
    % the relative
    pg_norm_g_rlt(i) = norm(pg_vec_g(:,i))/pg_norm(i);
    pg_H_norm_g_rlt(i) = norm(pg_vec_H_g(:,i))/pg_H_norm(i);
    % the relative
    pg_norm_c_rlt(i) = norm(pg_vec_c(:,i))/pg_norm(i);
    pg_H_norm_c_rlt(i) = norm(pg_vec_H_c(:,i))/pg_H_norm(i);
    % the relative
    pg_norm_h_rlt(i) = norm(pg_vec_h(:,i))/pg_norm(i);
    pg_H_norm_h_rlt(i) = norm(pg_vec_H_h(:,i))/pg_H_norm(i);
end
% look for the top 5 influenced edges for each edge absolutely
[pg_norm_maxk,ind_norm] = maxk(pg_norm,5);
[pg_H_norm_maxk,ind_H_norm] = maxk(pg_H_norm,5);
[pg_norm_g_maxk,ind_norm_g] = maxk(pg_norm_g,5);
[pg_H_norm_g_maxk,ind_H_norm_g] = maxk(pg_H_norm_g,5);
[pg_norm_c_maxk,ind_norm_c] = maxk(pg_norm_c,5);
[pg_H_norm_c_maxk,ind_H_norm_c] = maxk(pg_H_norm_c,5);
[pg_norm_h_maxk,ind_norm_h] = maxk(pg_norm_h,5);
[pg_H_norm_h_maxk,ind_H_norm_h] = maxk(pg_H_norm_h,5);
% look for the top 5 influenced edges for each edge relatively
[pg_norm_g_rlt_maxk,ind_norm_g_rlt] = maxk(pg_norm_g_rlt,5);
[pg_H_norm_g_rlt_maxk,ind_H_norm_g_rlt] = maxk(pg_H_norm_g_rlt,5);
[pg_norm_c_rlt_maxk,ind_norm_c_rlt] = maxk(pg_norm_c_rlt,5);
[pg_H_norm_c_rlt_maxk,ind_H_norm_c_rlt] = maxk(pg_H_norm_c_rlt,5);
[pg_norm_h_rlt_maxk,ind_norm_h_rlt] = maxk(pg_norm_h_rlt,5);
[pg_H_norm_h_rlt_maxk,ind_H_norm_h_rlt] = maxk(pg_H_norm_h_rlt,5);

% let us investigate the structure of these edges, if they contribute to a
% triangle

for i = 1:length(ind_norm)
    ind_norm_node1(i) = find(B1(:,ind_norm(i))==1);
    ind_norm_node2(i) = find(B1(:,ind_norm(i))==-1);
    ind_norm_tri(i) = nnz(B2t(:,ind_norm(i))) ;
end
for i = 1:length(ind_norm_g)
    ind_norm_g_node1(i) = find(B1(:,ind_norm_g(i))==1);
    ind_norm_g_node2(i) = find(B1(:,ind_norm_g(i))==-1);
    ind_norm_g_tri(i) = nnz(B2t(:,ind_norm_g(i))) ;
    ind_norm_g_rlt_tri(i) = nnz(B2t(:,ind_norm_g_rlt(i))) ;
end
for i = 1:length(ind_norm_c)
    ind_norm_c_node1(i) = find(B1(:,ind_norm_c(i))==1);
    ind_norm_c_node2(i) = find(B1(:,ind_norm_c(i))==-1);
    ind_norm_c_tri(i) = nnz(B2t(:,ind_norm_c(i))) ; 
    ind_norm_c_rlt_tri(i) = nnz(B2t(:,ind_norm_c_rlt(i))) ;
end
for i = 1:length(ind_norm_h)
    ind_norm_h_node1(i) = find(B1(:,ind_norm_h(i))==1);
    ind_norm_h_node2(i) = find(B1(:,ind_norm_h(i))==-1);
    ind_norm_h_tri(i) = nnz(B2t(:,ind_norm_h(i))) ;
    ind_norm_h_rlt_tri(i) = nnz(B2t(:,ind_norm_h_rlt(i))) ;
end

%% illustration
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.00], [0.08 0.02], [0.05 0.00]);
if ~make_it_tight,  clear subplot;  end
figure;
subplot(2,2,2);
color = string(['k','g','b','r']')';
ind_pg_norm = [ind_H_norm' ind_H_norm_g' ind_H_norm_c' ind_H_norm_h'];
ind_pg_norm_rlt = [ind_H_norm' ind_H_norm_g_rlt' ind_H_norm_c_rlt' ind_H_norm_h_rlt'];
label = "";
plt_network_page_rank(B1,ind_pg_norm,...
    repmat(color,length(ind_pg_norm),1),label,x_coord,y_coord);
axis tight
hold on;

subplot(2,2,4);
label = "";
plt_network_page_rank(B1,ind_pg_norm_rlt,...
    repmat(color,length(ind_pg_norm),1),label,x_coord,y_coord);
axis tight
%%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.08], [0.1 0.03], [0.08 0.01]);
% if ~make_it_tight,  clear subplot;  end
% figure;

subplot(2,2,1);
% for subplot 1
%scatter(ind_norm,pg_norm(ind_norm),'sk'); hold on;
scatter(ind_H_norm,pg_H_norm(ind_H_norm),100,'sk','filled'); hold on;
%scatter(ind_norm_g,pg_norm_g(ind_norm_g),'sg'); hold on;
scatter(ind_H_norm_g,pg_H_norm_g(ind_H_norm_g),100,'sg','filled'); hold on;
%scatter(ind_norm_c,pg_norm_c(ind_norm_c),'sb'); hold on;
scatter(ind_H_norm_c,pg_H_norm_c(ind_H_norm_c),100,'sb','filled'); hold on
%scatter(ind_norm_h,pg_norm_h(ind_norm_h),'sr'); hold on;
scatter(ind_H_norm_h,pg_H_norm_h(ind_H_norm_h),100,'sr','filled');
legend('Norm of PageRank vec', 'Norm of grad comp', 'Norm of curl comp', ...
    'Norm of harmonic comp');
xlim([1,130])
set(gca,'fontsize',12);
 
xticklabels({})
subplot(2,2,3);
% for subplot 2
%scatter(ind_norm,pg_norm(ind_norm)./pg_norm(ind_norm),'sk'); hold on;
scatter(ind_H_norm,pg_H_norm(ind_H_norm)./pg_H_norm(ind_H_norm),100,'sk','filled'); hold on;
scatter(ind_H_norm_g_rlt,pg_H_norm_g_rlt(ind_H_norm_g_rlt),100,'sg','filled'); hold on;
scatter(ind_H_norm_c_rlt,pg_H_norm_c_rlt(ind_H_norm_c_rlt),100,'sb','filled'); hold on
scatter(ind_H_norm_h_rlt,pg_H_norm_h_rlt(ind_H_norm_h_rlt),100,'sr','filled');
legend('rlt norm of PageRank vec', 'rlt norm of grad comp',...
    'rlt norm of curl comp', ...
    'rlt norm of harmonic comp');
xlim([1,130])
set(gca,'fontsize',12);
xlabel('Edge Index');

%% choose some edges, plot their page rank results
pos = [19,27,45,96];
% 19: large gradient
% 27: big harmonic
% 45: large influence on the rest
% 30: num_lower_nbrs large
% 96: num_upper_nbrs large, large curl component
 label = string(["19",'27','45','96']');
 %label=num2str(pos);
width_ratio = 10;
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.08]);
if ~make_it_tight,  clear subplot;  end
figure;
for i = 1:length(pos)
    subplot(2,2,i)
    h=plt_network(B1,pg_vec(:,pos(i)),x_coord,y_coord,width_ratio,1);
    axis tight
    node_label1(i) = find(B1(:,pos(i))==1);
    node_label2(i) = find(B1(:,pos(i))==-1);
    labeledge(h,node_label1(i),node_label2(i),label(i));
    h.EdgeFontSize=12;
end

