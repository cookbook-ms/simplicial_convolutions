% Script to create street network example analysis for London street
% networks
clc;
clear;
close all;
addpath('/Users/yangmaosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/London_street')
addpath('/Users/yangmaosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/simplicial-regularization')
%addpath '/tudelft.net/staff-bulk/ewi/insy/MMC/maosheng/simplicial-regularization/'
rng(123)
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
I = eye(num_edges);
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
% graph Laplacian
L0 = B1*B1';
% triangle Laplacian
L2 = B2'*B2;

% spectrum
[u1l, lam_l]  = eig(L1l); eig_L1l = diag(lam_l); u_g = u1l(:,50:end); lam_g = eig_L1l(50:end);
[u1, lam]= eig(L1); eig_L1 = diag(lam); u_h = u1(:,1:37);
[u1u, lam_u]= eig(L1u); eig_L1u = diag(lam_u); u_c = u1u(:,119:end);lam_c = eig_L1u(119:end);

% note that u_c don't necessarily have the same ordering as in u1;
%% Create initial flows and noisy realizations
% random harmonic initial condition
f_harm = P_h*randn(length(L1),1);
f_harm = f_harm/norm(f_harm);
% create synthetic low-gradient components
mu_g = 2;
[f_g,f_g_tilde] = create_low_grad_comp(lam_g,mu_g,u_g);
% f_g = (I+0.5*L1l)\randn(num_edges,1);
% f_g = f_g/norm(f_g);
% create synthetic low-curl components
mu_c = 5;
[f_c,f_c_tilde] = create_low_curl_comp(lam_c,mu_c,u_c);
% f_c = (I+0.5*L1u)\randn(num_edges,1);
% f_c = f_c/norm(f_c);

% noise free flow
f_noise_free = 10*f_harm ;%+ 5*f_g + 5*f_c  ;
q_noise_free= norm(B2'*f_noise_free);
p_noise_free= norm(B1*f_noise_free);

% noisy flow
% I think to show the ability of the new regularization framework we have
% to have the noise which includes different levels of the curl and
% gradient component, so there is inhomogenunity present in the noise
noise_g = B1'*randn(num_nodes,1); div_noise_g = norm(B1*noise_g);
noise_c = B2*randn(num_triangles,1); curl_noise_g = norm(B2'*noise_c);
noise = 0.25*noise_g + 0.25*noise_c;

f_initial = 1*f_noise_free ;%+ noise;
q_init = norm(B2'*f_initial); %curl
p_init = norm(B1*f_initial); %div

flows_grad = P_g*f_initial; div = B1*f_initial; div_total = norm(div);
flows_harm = P_h*f_initial;
flows_curl = P_c*f_initial; curl = B2'*f_initial; curl_total = norm(curl);

err_initial = norm(f_initial-f_noise_free)/norm(f_noise_free);
%% study the spectrum of the underlying edge flow
f_tilde_noise_free = u1'*f_noise_free; % f_noise_free;
f_tilde_noise = u1'*f_initial;

%%
ratio = 0.6;%0.2:0.05:1.0;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of sampled edges
    for ii = 1:1
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
        % the labeled
        f_in = f_initial.*mask;
        f_tilde_in = u1'*f_in;
        corr_in_single(i,ii) = corr(f_noise_free,f_in);

        % building the sampling matrix, which contains the same info as the
        % expanding matrix 
        sampling_matrix = zeros(M,num_edges);
        sampling_matrix_1 = diag(mask);
        sampling_matrix = sampling_matrix_1(any(sampling_matrix_1,2),:);
        mu = 0.5; alpha=0.17;
        % the MKL method
        [f_mkl,error] = mkl_interpolate(L1l,L1u,f_in,f_noise_free,mu,alpha,sampling_matrix_1);
        corr_out_single(i,ii) = corr(f_noise_free,f_mkl);
        f_tilde_mkl = u1'*f_mkl;

        % the approach in Jia2019
        mask_un = not(mask);
        % build the expanding matrix
        expanding_matrix = zeros(num_edges,nnz(mask_un));
        row_ind = find(mask_un==1);
        for j = 1:length(row_ind)
            expanding_matrix(row_ind(j),j) = 1;
        end
        % pseudo inverse method, which leads to bad performance, slightly
        % better than zero-fill
        f_unlabeled = pinv([B1*expanding_matrix;0.1*eye(nnz(mask_un))])...
            *[-B1*f_in;zeros(nnz(mask_un),1)];
        %f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
         %   [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_rgl_1 = f_in + expanding_matrix*f_unlabeled;
        f_tilde_rgl_1 = u1'*f_rgl_1;
        corr_rgl_1(i,ii) = corr(f_noise_free,f_rgl_1);
        
        % the approach in Schaub2021
        f_unlabeled_2 = pinv([B1*expanding_matrix;0.1*eye(nnz(mask_un));B2'*expanding_matrix])...
            *[-B1*f_in;zeros(nnz(mask_un),1);-B2'*f_in];
        f_rgl_2 = f_in + expanding_matrix*f_unlabeled_2;
        f_tilde_rgl_2 = u1'*f_rgl_2;
        corr_rgl_2(i,ii) = corr(f_noise_free,f_rgl_2);
     % our approach in jia2019 fashion
        f_unlabeled_2 = pinv([0.6*B1*expanding_matrix;0.1*eye(nnz(mask_un));0.3*B2'*expanding_matrix])...
            *[0.6*-B1*f_in;zeros(nnz(mask_un),1);0.3*-B2'*f_in];
        f_mkl_2 = f_in + expanding_matrix*f_unlabeled_2;
        corr_mkl_2(i,ii) = corr(f_noise_free,f_mkl_2);
        %% consider l1 denoising
        cvx_begin
        variables f_opt(num_edges)  ;
        minimize(0.0* norm(f_opt)+0.5*norm(B1*f_opt,1)+0.5*norm(B2'*f_opt,1));
        subject to
                sampling_matrix*f_in == sampling_matrix*f_opt;
        
        cvx_end
        f_l1 = f_opt;
        f_tilde_l1 = u1'*f_l1;
        corr_l1(i,ii) = corr(f_noise_free,f_l1);
        err_l1(i,ii) = norm(f_opt-f_noise_free)/norm(f_noise_free);
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out_mean(i) = mean(corr_out_single(i,:),2);
    corr_rgl_1_mean(i) = mean(corr_rgl_1(i,:),2);
    corr_rgl_2_mean(i) = mean(corr_rgl_2(i,:),2);
    corr_mkl_2_mean(i) = mean(corr_mkl_2(i,:),2);
    corr_l1_mean(i) = mean(corr_l1(i,:),2);
    
end

%%
figure;
plot(ratio,corr_in); hold on;
plot(ratio,corr_out_mean); hold on;
plot(ratio,corr_rgl_1_mean); hold on;
plot(ratio,corr_rgl_2_mean); hold on;
plot(ratio,corr_mkl_2_mean); hold on;
plot(ratio,corr_l1_mean)
legend('Initial','MKL','RGL 1', 'RGL 2', 'MKL 2','l1')
