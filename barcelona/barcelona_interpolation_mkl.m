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
t4 = readtable('flow_vec.csv');
num_nodes = 1020; num_edges = 1798; num_triangles = 164;
I = eye(num_edges);

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
% note that u_c don't necessarily have the same ordering as in u1;
%% Create initial flows and noisy realizations
f_noise_free = f;%/norm(f);
q_noise_free= norm(B2'*f_noise_free); 
p_noise_free= norm(B1*f_noise_free);
noise_g = B1'*randn(num_nodes,1); div_noise_g = norm(B1*noise_g);
noise_c = B2*randn(num_triangles,1); curl_noise_g = norm(B2'*noise_c);
noise = 0.25*noise_g + 0.25*noise_c;

% noisy flow
f_initial = 1*f_noise_free ;%+ 0.05*noise;
q_init = norm(B2'*f_initial); %curl
p_init = norm(B1*f_initial); %div

flows_grad = P_g*f_initial; div = B1*f_initial; div_total = norm(div);
flows_harm = P_h*f_initial;
flows_curl = P_c*f_initial; curl = B2'*f_initial; curl_total = norm(curl);
err_initial = norm(f_initial-f_noise_free)/norm(f_noise_free);
%% study the spectrum of the underlying edge flow
% f_tilde_noise_free = u1'*f_noise_free; % f_noise_free;
% f_tilde_noise = u1'*f_initial;
% figure;
% plot(eig_L1,f_tilde_noise_free);
% hold on ;
% plot(eig_L1, f_tilde_noise);
%%
ratio = 0.2:0.1:1.0;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of sampled edges
    for ii = 1:1
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
        % the labeled
        f_in = f_initial.*mask;
        corr_in_single(i,ii) = corr(f_noise_free,f_in);

        % building the sampling matrix, which contains the same info as the
        % expanding matrix 
        sampling_matrix = zeros(M,num_edges);
        sampling_matrix_1 = diag(mask);
        sampling_matrix = sampling_matrix_1(any(sampling_matrix_1,2),:);
        alpha=0.35; mu = 0.3; 
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
        f_unlabeled = pinv([B1*expanding_matrix;0.1*eye(nnz(mask_un))])...
            *[-B1*f_in;zeros(nnz(mask_un),1)];
       % f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
        %    [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_rgl_1 = f_in + expanding_matrix*f_unlabeled;
        corr_rgl_1(i,ii) = corr(f_noise_free,f_rgl_1);
        
        % the approach in Schaub2021
        f_unlabeled_2 = pinv([B1*expanding_matrix;0.1*eye(nnz(mask_un));B2'*expanding_matrix])...
            *[-B1*f_in;zeros(nnz(mask_un),1);-B2'*f_in];
        f_rgl_2 = f_in + expanding_matrix*f_unlabeled_2;
        corr_rgl_2(i,ii) = corr(f_noise_free,f_rgl_2);

        % our approach in jia2019 fashion
        f_unlabeled_2 = pinv([0.8*B1*expanding_matrix;0.4*eye(nnz(mask_un));0.6*B2'*expanding_matrix])...
            *[0.8*-B1*f_in;zeros(nnz(mask_un),1);0.6*-B2'*f_in];
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
