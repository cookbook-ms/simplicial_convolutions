% Script to create street network example analysis for London street
% networks
clc;
clear;
close all;
addpath('/Users/maosheng/OneDrive - Delft University of Technology/Topological filter/codes from 15 update/experiments/simplicial-regularization')
%addpath '/tudelft.net/staff-bulk/ewi/insy/MMC/maosheng/simplicial-regularization/'
%rng(1223)
t1 = readtable('B1.csv');
t2 = readtable('B2t.csv');
t3 = readtable('coordinate.csv');
t4 = readtable('flow_vec.csv');
num_nodes = 546; num_edges = 1088; num_tri = 112;
I = eye(num_edges);
coordinate = t3{:,:}';
B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t';
L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;
f = t4{:,:};
f = f - mean(abs(f));
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
[u1l, lam_l]  = eig(L1l); eig_L1l = diag(lam_l); u_g = u1l(:,544:end); lam_g = eig_L1l(544:end);
[u1, lam]= eig(L1); eig_L1 = diag(lam); u_h = u1(:,1:431);
[u1u, lam_u]= eig(L1u); eig_L1u = diag(lam_u); u_c = u1u(:,977:end);lam_c = eig_L1u(977:end);

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
f_noise_free = 1*f_harm + 5*f_g + 5*f_c ;
%
f_noise_free = f;%/norm(f);
q_noise_free= norm(B2'*f_noise_free);
p_noise_free= norm(B1*f_noise_free);
% noisy flow
f_initial = 1*f_noise_free + 1000*randn(num_edges,1);
q_init = norm(B2'*f_initial); %curl
p_init = norm(B1*f_initial); %div

flows_grad = P_g*f_initial; div = B1*f_initial; div_total = norm(div);
flows_harm = P_h*f_initial;
flows_curl = P_c*f_initial; curl = B2'*f_initial; curl_total = norm(curl);
err_initial = norm(f_initial-f_noise_free)/norm(f_noise_free);
%% study the spectrum of the underlying edge flow
f_tilde_noise_free = u1'*f_noise_free; % f_noise_free;
f_tilde_noise = u1'*f_initial;
figure;
plot(eig_L1,f_tilde_noise_free)

%%
ratio = 0.2:0.1:1.0;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of sampled edges
    mask = zeros(num_edges,1);
    mask(randperm(numel(mask), M)) = 1;
    for ii = 1:10
        % the labeled
        f_in = f_initial.*mask;
        corr_in_single(i,ii) = corr(f_noise_free,f_in);

        % building the sampling matrix, which contains the same info as the
        % expanding matrix
        sampling_matrix = zeros(M,num_edges);
        sampling_matrix_1 = diag(mask);
        sampling_matrix = sampling_matrix_1(any(sampling_matrix_1,2),:);
        mu = 0.2; alpha=0.81;
        % the MKL method
        [f_mkl,error] = mkl_interpolate(L1l,L1u,f_in,f_noise_free,mu,alpha,sampling_matrix_1);
        corr_out_single(i,ii) = corr(f_noise_free,f_mkl);

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
        f_unlabeled = pinv([B1*expanding_matrix;0.5*eye(nnz(mask_un))])...
            *[-B1*f_in;zeros(nnz(mask_un),1)];
        %f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
        %   [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_rgl_1 = f_in + expanding_matrix*f_unlabeled;
        corr_rgl_1(i,ii) = corr(f_noise_free,f_rgl_1);

        % the approach in Schaub2021
        f_unlabeled_2 = pinv([B1*expanding_matrix;0.3*eye(nnz(mask_un));B2'*expanding_matrix])...
            *[-B1*f_in;zeros(nnz(mask_un),1);-B2'*f_in];
        f_rgl_2 = f_in + expanding_matrix*f_unlabeled_2;
        corr_rgl_2(i,ii) = corr(f_noise_free,f_rgl_2);
            % our approach in jia2019 fashion
        f_unlabeled_2 = pinv([0.6*B1*expanding_matrix;0.3*eye(nnz(mask_un));0.3*B2'*expanding_matrix])...
            *[0.6*-B1*f_in;zeros(nnz(mask_un),1);0.3*-B2'*f_in];
        f_mkl_2 = f_in + expanding_matrix*f_unlabeled_2;
        corr_mkl_2(i,ii) = corr(f_noise_free,f_mkl_2);
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out_mean(i) = mean(corr_out_single(i,:),2);
    corr_rgl_1_mean(i) = mean(corr_rgl_1(i,:),2);
    corr_rgl_2_mean(i) = mean(corr_rgl_2(i,:),2);
    corr_mkl_2_mean(i) = mean(corr_mkl_2(i,:),2);
end

%%
figure;
plot(ratio,corr_in); hold on;
plot(ratio,corr_out_mean); hold on;
plot(ratio,corr_rgl_1_mean); hold on; 
plot(ratio,corr_rgl_2_mean); hold on;
plot(ratio,corr_mkl_2_mean);
legend('Initial','MKL','RGL 1', 'RGL 2', 'MKL 2')
