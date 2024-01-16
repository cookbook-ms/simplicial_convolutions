clc;
clear;
close all;
addpath('/Users/maosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/forex')
%rng(1223)
t1 = readtable('B1_FX_1538755200.csv');
t2 = readtable('B2t_FX_1538755200.csv');
t4 = readtable('flow_FX_1538755200.csv');
B1 = t1{:,:};
B2t = t2{:,:};
B2 = B2t'; 
% the edge flow is the mid of the ask and bid prices
f = t4{:,1};
% Hodge Laplacian
L1l = B1'*B1; L1u = B2*B2';
L1 = L1l + L1u;
% edge Laplacian
Le = L1l;
% graph Laplacian
L0 = B1*B1';
% triangle Laplacian
L2 = B2'*B2;
num_nodes = size(B1,1); num_edges = size(B1,2); num_tri = size(B2,2);
I = eye(num_edges);
[u1, lam]= eig(L1); eig_L1 = diag(lam);
% check the divergence and curl of the original signal
div_f = norm(B1*f); % the actual transition signal is div-free
curl_f = norm(B2'*f); % it is not curl-free
%% interpolation task
ratio = 0.05:0.05:1;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    for ii = 1:2
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
        sampling_matrix = zeros(M,num_edges);
        sampling_matrix_1 = diag(mask);
        sampling_matrix = sampling_matrix_1(any(sampling_matrix_1,2),:);
        % the labeled
        f_in = f.*mask;
        corr_in_single(i,ii) = corr(f,f_in);
        % the filter method
        f_filtered1 = I*f_in;
        corr_out_single(i,ii) = corr(f,f_filtered1);
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
        f_unlabeled = pinv([B2'*expanding_matrix;2*eye(nnz(mask_un))])...
            *[-B2'*f_in;zeros(nnz(mask_un),1)];
        %f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
        %    [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_jia2019 = f_in + expanding_matrix*f_unlabeled;
        corr_jia2019(i,ii) = corr(f,f_jia2019);
        curl_jia2019(i,ii) = norm(B2'*f_jia2019);
        count_jia2019(i,ii) =  nnz(abs(f-f_jia2019) < 0.1*f)/num_edges;
        %% l1 norm based regularization 
        % here we shall the constraint to preserve the observed labels
        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 2* norm(B2'*L1u^0*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_1 = f_opt;
        err_l1_1 = norm(f_l1_1-f)/norm(f);
        corr_l1_1(i,ii) = corr(f,f_l1_1);
        % consider other metrics 
        count_l1_1(i,ii) = nnz(abs(f-f_l1_1) < 0.1*f)/num_edges;
        curl_l1_1(i,ii) = norm(B2'*f_l1_1); 

        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 2* norm(B2'*L1u^1*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_2 = f_opt;
        err_l1_2 = norm(f_l1_2-f)/norm(f);
        corr_l1_2(i,ii) = corr(f,f_l1_2);
        % consider other metrics 
        count_l1_2(i,ii) = nnz(abs(f-f_l1_2) < 0.1*f)/num_edges;
        curl_l1_2(i,ii) = norm(B2'*f_l1_2); 

        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 2* norm(L1u*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_3 = f_opt;
        err_l1_3 = norm(f_l1_3-f)/norm(f);
        corr_l1_3(i,ii) = corr(f,f_l1_3);
        % consider other metrics 
        count_l1_3(i,ii) = nnz(abs(f-f_l1_3) < 0.1*f)/num_edges;
        curl_l1_3(i,ii) = norm(B2'*f_l1_3); 
        
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out(i) = mean(corr_out_single(i,:),2);
    corr_jia2019_mean(i) = mean(corr_jia2019(i,:),2);
    corr_l1_1_mean(i) = mean(corr_l1_1(i,:),2);
    corr_l1_2_mean(i) = mean(corr_l1_2(i,:),2);
    corr_l1_3_mean(i) = mean(corr_l1_3(i,:),2);

    count_l1_1_mean(i) = mean(count_l1_1(i,:),2);
    count_l1_2_mean(i) = mean(count_l1_2(i,:),2);
    count_l1_3_mean(i) = mean(count_l1_3(i,:),2);
    count_jia2019_mean(i) = mean(count_jia2019(i,:),2);

    curl_l1_1_mean(i) = mean(curl_l1_1(i,:),2);
    curl_l1_2_mean(i) = mean(curl_l1_2(i,:),2);
    curl_l1_3_mean(i) = mean(curl_l1_3(i,:),2);
    curl_jia2019_mean(i) = mean(curl_jia2019(i,:),2);
end
%%
% filename = 'forex_interpolation';
% save(filename)

subplot(2,2,3);
plot(ratio,corr_jia2019_mean,'--','LineWidth',3.5); hold on;
plot(ratio,corr_l1_mean,'LineWidth',2); hold on;
% plot(ratio,corr_l1_2_mean,'LineWidth',2); hold on;
% plot(ratio,corr_l1_3_mean,'LineWidth',2); hold on;

plot(ratio,corr_in,'k','LineWidth',2); hold on;
legend( 'ssl',  'STF-0','STF-2','STF-1','zero fill')
set(gca,'fontsize',14)
xlabel('Ratio labeled')
grid on;

subplot(2,2,4); 
plot(ratio, curl_jia2019_mean,'--','LineWidth',2.5); hold on; 
plot(ratio, curl_l1_mean,'LineWidth',2);
% plot(ratio, curl_l1_2_mean,'LineWidth',2);
% plot(ratio, curl_l1_3_mean,'LineWidth',2);
legend( 'ssl',  'STF-0','STF-2','STF-1');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('Ratio labeled')
grid on;
%% plot the figures of the results from the hpc
addpath('/Users/maosheng/Documents/subtightplot')
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.08], [0.1 0.03], [0.08 0.01]);
if ~make_it_tight,  clear subplot;  end
 
subplot(2,2,3);
plot(ratio,corr_jia2019,'--','LineWidth',3.5); hold on;
plot(ratio,corr_l1_1_mean,'LineWidth',2); hold on;
plot(ratio,corr_l1_4_mean,'LineWidth',2); hold on;
plot(ratio,corr_l1_2_mean,'LineWidth',2); hold on;
plot(ratio,corr_l1_5_mean,'LineWidth',2); hold on;
plot(ratio,corr_l1_3_mean,'LineWidth',2); hold on;
plot(ratio,corr_out,'k','LineWidth',2)

legend('$\ell_2$ based', '0-STF','1-STF',...
    '2-STF','3-STF','4-STF',....
    'zero fill','Interpreter','latex');

set(gca,'fontsize',14)
xlabel('SNR')

subplot(2,2,4);
plot(ratio,curl_jia2019_mean,'--','LineWidth',2); hold on;
plot(ratio, curl_l1_1_mean,'LineWidth',2)
plot(ratio, curl_l1_4_mean,'LineWidth',2)
plot(ratio, curl_l1_2_mean,'LineWidth',2)
plot(ratio, curl_l1_5_mean,'LineWidth',2)
plot(ratio, curl_l1_3_mean,'LineWidth',2)

legend('$\ell_2$ based', '0-STF','1-STF',...
    '2-STF','3-STF','4-STF','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')