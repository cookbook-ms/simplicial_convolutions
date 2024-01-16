clear all
% this is a script to process flows on Chicago transportation network
addpath('~/Documents/chebfun')
addpath('/Users/maosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/Lastfm-dataset-1k')
cd /Users/maosheng/Documents/cvx
% cvx_setup

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
B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t'; B1*B2;
L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;
I = eye(num_edges);
% the original data can be used for interpolation, where the filtering
% method does not perform well though
f = t3{:,:};
% check the divergence and curl of the original signal
div_f = norm(B1*f); % the actual transition signal is div-free
curl_f = norm(B2'*f); % it is not curl-free

%% generate noisy flow
% rng(12);
% noise = randn(num_edges,1);
% f_noisy = f + noise;
% div_noisy = norm(B1*f_noisy);
% curl_noisy = norm(B2'*f_noisy);
% err_noisy = norm(f_noisy-f)/norm(f);
% %% perform l2 denoising
% mu = 0.5;
% f_l2 = (I+mu*L1l)\f_noisy;
% err_l2 = norm(f_l2-f)/norm(f);
% div_l2 = norm(B1*f_l2);
% curl_l2 = norm(B2'*f_l2);
% %% consider l1 denoising
% cvx_begin
%     variables f_opt(num_edges) z(num_nodes);
%     minimize( norm(f_noisy-f_opt)+0.08*norm(z,1) );
%     subject to 
%         B1*f_opt == z;
% cvx_end
% f_l1 = f_opt;
% err_l1 = norm(f_opt-f)/norm(f);
% div_l1 = norm(B1*f_l1);
% curl_l1 = norm(B2'*f_l1);

%%
num_realizations = 1;
snr_db = -12:1:-11;
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
        % grid search mu
%         mu = 0.1:0.1:1; %0.3;
%         for iii = 1:length(mu)
%             f_l2 = (I+mu(iii)*L1l)\f_noisy;
%             err_l2(i,ii,iii) = norm(f_l2-f)/norm(f);
%         end
%         [~,M(i,ii)] = min(err_l2(i,ii,:));

        mu = 0.5; % this is the optimal mu when snr is 3db
        f_l2 = (I+mu*L1l)\f_noisy;
        err_l2(i,ii) = norm(f_l2-f)/norm(f);
        div_l2(i,ii) = norm(B1*f_l2);
        curl_l2(i,ii) = norm(B2'*f_l2);
        %% consider l1 denoising
        cvx_begin
        variables f_opt(num_edges);
        minimize(1* norm(f_noisy-f_opt)+0.5*norm(B1*f_opt,1) );
        cvx_end
        f_l1_1 = f_opt;
        err_l1_1(i,ii) = norm(f_l1_1-f)/norm(f);
        div_l1_1(i,ii) = norm(B1*f_l1_1);
        curl_l1_1(i,ii) = norm(B2'*f_l1_1);

%         cvx_begin
%         variables f_opt(num_edges);
%         minimize(1* norm(f_noisy-f_opt)+0.5*norm(B1*L1l^2*f_opt,1) );
%         cvx_end
        f_l1_2 = f_opt;
        err_l1_2(i,ii) = norm(f_l1_2-f)/norm(f);
        div_l1_2(i,ii) = norm(B1*f_l1_2);
        curl_l1_2(i,ii) = norm(B2'*f_l1_2);

%         cvx_begin
%         variables f_opt(num_edges);
%         minimize(1* norm(f_noisy-f_opt)+0.5*norm(B1*L1l*f_opt,1) );
%         cvx_end
        f_l1_3 = f_opt;
        err_l1_3(i,ii) = norm(f_l1_3-f)/norm(f);
        div_l1_3(i,ii) = norm(B1*f_l1_3);
        curl_l1_3(i,ii) = norm(B2'*f_l1_3);
    end
end
%%
div_l1_1_mean = mean(div_l1_1,2);
div_l1_2_mean = mean(div_l1_2,2);
div_l1_3_mean = mean(div_l1_3,2);
div_l2_mean = mean(div_l2,2);

err_l1_1_mean = mean(err_l1_1,2);
err_l1_2_mean = mean(err_l1_2,2);
err_l1_3_mean = mean(err_l1_3,2);
err_l2_mean = mean(err_l2,2);
err_noisy_mean = mean(err_noisy,2);

% filename = 'lastfm_denoise';
% save(filename)
%%
figure;
subplot(2,2,1);
plot(snr_db,err_l2_mean,'--','LineWidth',3.5); hold on;
plot(snr_db,err_l1_mean,'LineWidth',2); hold on;
% plot(snr_db,err_l1_2_mean,'LineWidth',2); hold on;
% plot(snr_db,err_l1_3_mean,'LineWidth',2); hold on;
plot(snr_db,err_noisy_mean,'k','LineWidth',2)
xlim([-12,6])
legend('$\ell_2$ regularizer', 'STF-0','STF-2','STF-1','noisy','Interpreter','latex');
%set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')
grid on

subplot(2,2,2);
plot(snr_db,div_l2_mean,'LineWidth',2); hold on;
plot(snr_db, div_l1_mean,'LineWidth',2); hold on;
% plot(snr_db, div_l1_2_mean,'LineWidth',2); hold on;
% plot(snr_db, div_l1_3_mean,'LineWidth',2)
xlim([-12,6])
legend('$\ell_2$ regularizer', 'STF-0','STF-2','STF-1','Interpreter','latex');
set(gca, 'YScale', 'log')
set(gca,'fontsize',14)
xlabel('SNR')
grid on

