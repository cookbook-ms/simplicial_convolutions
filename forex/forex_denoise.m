clc;
clear;
close all;
addpath('/Users/mmstudio/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/forex')
%rng(1223)
t1 = readtable('B1_FX_1538755200.csv');
t2 = readtable('B2t_FX_1538755200.csv');
t4 = readtable('flow_FX_1538755200.csv');
B1 = t1{:,:};
B2t = t2{:,:};
B2 = B2t';
num_nodes = size(B1,1); num_edges = size(B1,2); num_tri = size(B2,2);
I = eye(num_edges);
% the edge flow is the mid of the ask and bid prices
f = t4{:,2};
% Hodge Laplacian
L1l = B1'*B1; L1u = B2*B2';
L1 = L1l + L1u;
% edge Laplacian
Le = L1l;
% graph Laplacian
i0 = I;
i0(5,5)=2;
L0 = B1*i0*B1';
[u0, lam_0]= eig(L0);
% triangle Laplacian
L2 = B2'*B2;

[u1, lam]= eig(L1); eig_L1 = diag(lam);
% check the divergence and curl of the original signal
div_f = norm(B1*f); % the actual transition signal is div-free
curl_f = norm(B2'*f); % it is not curl-free
f_proj = 1/num_nodes*L1l*f;
%%
num_realizations = 1;
snr_db = -3 %-6:1:12;
snr = 10.^(snr_db/10);
power_flow = norm(f,2);
power_noise = power_flow./snr/num_edges;
for i = 1:length(snr)
    for ii = 1:num_realizations
        noise = power_noise(i)*randn(num_edges,1);
        f_noisy = f + noise;
        div_noisy(i,ii) = norm(B1*f_noisy);
        curl_noisy(i,ii) = norm(B2'*f_noisy);
        err_noisy(i,ii) = norm(f_noisy-f)/norm(f);

        %% perform l2 denoising
%         mu = 0.05:0.1:5; %0.3;
%         for iii = 1:length(mu)
%             f_l2(:,iii) = (I+mu(iii)*L1u)\f_noisy;
%             err_l2(i,ii,iii) = norm(f_l2(:,iii)-f)/norm(f);
%         end
%         [~,M(i,ii)] = min(err_l2(i,ii,:));
        mu = 2; % this is the optimal mu when snr is 3db 
        f_l2 = (I+mu*L1u)\f_noisy;
        err_l2(i,ii) = norm(f_l2-f)/norm(f);
        div_l2(i,ii) = norm(B1*f_l2);
        curl_l2(i,ii) = norm(B2'*f_l2);
        %% consider l1 denoising
        cvx_begin
        variables f_opt(num_edges);
        minimize(1* norm(f_noisy-f_opt)+2*norm(B2'*f_opt,1) );
        cvx_end
        f_l1_1 = f_opt;
        err_l1_1(i,ii) = norm(f_l1_1-f)/norm(f);
        div_l1_1(i,ii) = norm(B1*f_l1_1);
        curl_l1_1(i,ii) = norm(B2'*f_l1_1);

        cvx_begin
        variables f_opt(num_edges);
        minimize(1* norm(f_noisy-f_opt)+2*norm(B2'*L1u^2*f_opt,1) );
        cvx_end
        f_l1_2 = f_opt;
        err_l1_2(i,ii) = norm(f_l1_2-f)/norm(f);
        div_l1_2(i,ii) = norm(B1*f_l1_2);
        curl_l1_2(i,ii) = norm(B2'*f_l1_2);

        cvx_begin
        variables f_opt(num_edges);
        minimize(1* norm(f_noisy-f_opt)+2*norm(L1u^2*f_opt,1) );
        cvx_end
        f_l1_3 = f_opt;
        err_l1_3(i,ii) = norm(f_l1_3-f)/norm(f);
        div_l1_3(i,ii) = norm(B1*f_l1_3);
        curl_l1_3(i,ii) = norm(B2'*f_l1_3);
    end
end
%%
curl_l1_1_mean = mean(curl_l1_1,2);
curl_l1_2_mean = mean(curl_l1_2,2);
curl_l1_3_mean = mean(curl_l1_3,2);
curl_l2_mean = mean(curl_l2,2);

err_l1_1_mean = mean(err_l1_1,2);
err_l1_2_mean = mean(err_l1_2,2);
err_l1_3_mean = mean(err_l1_3,2);
err_l2_mean = mean(err_l2,2);
err_noisy_mean = mean(err_noisy,2);

%%
figure;
subplot(2,2,1);
plot(snr_db,err_l2_mean,'--','LineWidth',3.5); hold on;
plot(snr_db,err_l1_mean,'LineWidth',2); hold on;
% plot(snr_db,err_l1_2_mean,'LineWidth',2); hold on;
% plot(snr_db,err_l1_3_mean,'LineWidth',2); hold on;
plot(snr_db,err_noisy_mean,'k','LineWidth',2)
xlim([-6,12])
legend('$\ell_2$ regularizer', 'STF-0','STF-2','STF-1','noisy','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')
grid on;

subplot(2,2,2);
plot(snr_db,curl_l2_mean,'LineWidth',2); hold on;
plot(snr_db, curl_l1_mean,'LineWidth',2); hold on
% plot(snr_db, curl_l1_2_mean,'LineWidth',2); hold on
% plot(snr_db, curl_l1_3_mean,'LineWidth',2); hold on
xlim([-6,12])
legend('$\ell_2$ regularizer', 'STF-0','STF-2','STF-1','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')
grid on;

%% plot the figures of the results from hpc
addpath('/Users/maosheng/Documents/subtightplot')
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.08], [0.1 0.03], [0.08 0.01]);
if ~make_it_tight,  clear subplot;  end
figure;
subplot(2,2,1);
plot(snr_db,err_l2_mean,'--','LineWidth',3.5); hold on;
plot(snr_db,err_l1_1_mean,'LineWidth',2); hold on;
plot(snr_db,err_l1_4_mean,'LineWidth',2); hold on;
plot(snr_db,err_l1_2_mean,'LineWidth',2); hold on;
plot(snr_db,err_l1_5_mean,'LineWidth',2); hold on;
plot(snr_db,err_l1_3_mean,'LineWidth',2); hold on;
plot(snr_db,err_noisy_mean,'k','LineWidth',2)
xlim([-6,12])
legend('$\ell_2$ based', '0-STF','1-STF',...
    '2-STF','3-STF','4-STF',....
    'noisy','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')

subplot(2,2,2);
plot(snr_db,curl_l2_mean,'--','LineWidth',2); hold on;
plot(snr_db, curl_l1_1_mean,'LineWidth',2)
plot(snr_db, curl_l1_4_mean,'LineWidth',2)
plot(snr_db, curl_l1_2_mean,'LineWidth',2)
plot(snr_db, curl_l1_5_mean,'LineWidth',2)
plot(snr_db, curl_l1_3_mean,'LineWidth',2)
xlim([-6,12])
legend('$\ell_2$ based', '0-STF','1-STF',...
    '2-STF','3-STF','4-STF','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')