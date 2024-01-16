clear all
% this is a script to process flows on Chicago transportation network
addpath('~/Documents/chebfun')
addpath('/Users/maosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/Lastfm-dataset-1k')
cd /Users/maosheng/Documents/cvx
cvx_setup
%% transition by song
t1 = readtable('B1.csv');
t2 = readtable('B2t.csv');
t3 = readtable('flow_vec.csv');
num_nodes = 3092; num_edges = 5507; num_tri = 571;

%% transition by artist
t1 = readtable('B1-artist.csv');
t2 = readtable('B2t-artist.csv');
t3 = readtable('flow_vec-artist.csv');
num_nodes = 657; num_edges = 1997; num_tri = 1276;

%% build in matlab
B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t';  
L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;
I = eye(num_edges);
% the original data can be used for interpolation, where the filtering
% method does not perform well though
f = t3{:,:};
% check the divergence and curl of the original signal
div_f = norm(B1*f); % the actual transition signal is div-free
curl_f = norm(B2'*f); % it is not curl-free

%% interpolation task
ratio = 0.05:0.05:1;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    for ii = 1:10
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
        f_unlabeled = pinv([B1*expanding_matrix;0.5*eye(nnz(mask_un))])...
            *[-B1*f_in;zeros(nnz(mask_un),1)];
        %f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
        %    [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_jia2019 = f_in + expanding_matrix*f_unlabeled;
        corr_jia2019(i,ii) = corr(f,f_jia2019);
        div_jia2019(i,ii) = norm(B1*f_jia2019);
        %% l1 norm based regularization 
        % here we shall the constraint to preserve the observed labels
        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 0.5* norm(B1*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_1 = f_opt;
        err_l1_1 = norm(f_l1_1-f)/norm(f);
        corr_l1_1(i,ii) = corr(f,f_l1_1);
        div_l1_1(i,ii) = norm(B1*f_l1_1); 

        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 0.5* norm(B1*L1l*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_2 = f_opt;
        err_l1_2 = norm(f_l1_2-f)/norm(f);
        corr_l1_2(i,ii) = corr(f,f_l1_2);
        div_l1_2(i,ii) = norm(B1*f_l1_2); 

        cvx_begin
            variable f_opt(num_edges);
            minimize( 1*norm(f_opt-f_in) + 0.5* norm(L1l*f_opt,1) );
            subject to
                sampling_matrix*f_opt == sampling_matrix*f_in;
        cvx_end
        f_l1_3 = f_opt;
        err_l1_3 = norm(f_l1_3-f)/norm(f);
        corr_l1_3(i,ii) = corr(f,f_l1_3);
        div_l1_3(i,ii) = norm(B1*f_l1_3); 
        
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out(i) = mean(corr_out_single(i,:),2);
    corr_jia2019_mean(i) = mean(corr_jia2019(i,:),2);
    corr_l1_1_mean(i) = mean(corr_l1_1(i,:),2);
    corr_l1_2_mean(i) = mean(corr_l1_2(i,:),2);
    corr_l1_3_mean(i) = mean(corr_l1_3(i,:),2);

    div_l1_1_mean(i) = mean(div_l1_1(i,:),2);
    div_l1_2_mean(i) = mean(div_l1_2(i,:),2);
    div_l1_3_mean(i) = mean(div_l1_3(i,:),2);
    div_jia2019_mean(i) = mean(div_jia2019(i,:),2);
end
%%
% filename = 'lastfm_interpolation';
% save(filename)

% %%
% figure;
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
plot(ratio, div_jia2019_mean,'--','LineWidth',2.5); hold on; 
plot(ratio, div_l1_mean,'LineWidth',2); hold on;
% plot(ratio, div_l1_2_mean,'LineWidth',2); hold on;
% plot(ratio, div_l1_3_mean,'LineWidth',2); hold on;
legend( 'ssl',  'STF-0','STF-2','STF-1');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('Ratio labeled')
grid on