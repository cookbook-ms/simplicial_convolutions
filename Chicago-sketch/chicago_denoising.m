% Script to create street network example analysis for London street
% networks
clc;
clear;
close all;
addpath('/Users/yangmaosheng/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Topological filter/codes from 15 update/experiments/Chicago-sketch')
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
plot(eig_L1,f_tilde_noise_free);
hold on ;
plot(eig_L1, f_tilde_noise);
%% solution comparison
% f_1 multi-kernel solution grid search: two parameters
% f_5 multi-kernel solution optimal parameter design
%
% f_2 edge laplacian solution grid search: one parameter
% f_3 hodge laplacian solution grid search: one parameter
% f_4 line graph solution grid search: one parameter
%%
%alpha = 0.001:0.01:0.99; mu = 1:1:100;
best = 1e6;
alpha = 0.01:0.01:1;
for i = 1:length(alpha)
    %alpha(i)=0.9322
    if alpha(i)==1
        err_1(i)=1;
        continue;
    end
    for mu =  0.1:0.1:1
        %mu = 0.1
        [f_1, err_1(i)] = mkl_denoise(L1l,L1u,f_initial,f_noise_free,mu,alpha(i));
        if err_1(i) < best
            best = err_1(i);
            alpha_best_grid = alpha(i);
            mu_best_grid = mu;
        end
    end
end
figure;
plot( 0.01:0.01:1,err_1);
err_1_best=best;
%f_1 = (I+ mu_best_grid/alpha_best_grid *L1l + mu_best_grid/(1-alpha_best_grid)*L1u)\f_initial;
f_1 = mkl_denoise(L1l,L1u,f_initial,f_noise_free,mu_best_grid,alpha_best_grid);
norm(f_1-f_noise_free)/norm(f_noise_free)
p_1= norm(B1*f_1); q_1 = norm(B2'*f_1); 
%%
best = 1e6;
rho_2 = 0.01:0.01:1;
for i = 1:length(rho_2)
    f_2 = (I+ rho_2 (i)*L1l )\f_initial;
    err_2(i) = norm(f_2-f_noise_free)/norm(f_noise_free);
    if  err_2(i) < best
        best =  err_2(i);
        rho_2_best = rho_2(i);
    end
end
err_2_best=best;
figure;
plot(rho_2,err_2);
f_2 = (I+ rho_2_best *L1l )\f_initial;
%% consider l1 denoising
cvx_begin
    variables f_opt(num_edges)  ;
    minimize( norm(f_initial-f_opt)+0.05*norm(B1*f_opt,1)+0.1*norm(B2'*f_opt,1));
 
cvx_end
f_l1 = f_opt;
err_l1 = norm(f_opt-f_noise_free)/norm(f_noise_free);
div_l1 = norm(B1*f_l1);
curl_l1 = norm(B2'*f_l1);

%%
best = 1e6;
rho_3=0.01:0.01:1;
for i = 1:length(rho_3)
    f_3 = (I+ rho_3(i) *L1l + rho_3(i) *L1u)\f_initial;
    err_3(i) = norm(f_3-f_noise_free)/norm(f_noise_free);
    if err_3(i) < best
        best = err_3(i);
        rho_3_best = rho_3(i) ;
    end
end
err_3_best=best;
figure;
f_3 = (I+ rho_3_best *L1l + rho_3_best*L1u)\f_initial;
plot(rho_3,err_3);
%% MKL - optimal parameter design
best = 1e6;
for mu_opt = 0.1:0.1:1%mu_best_grid
    mu = 0.2;
    alpha_opt = p_init/(p_init+q_init);
    alpha_opt = p_noise_free/(p_noise_free+q_noise_free)
    for t = 1:30 % number of iterations for alternating minimization
        [f_5,err_5(t)]= mkl_denoise(L1l,L1u,f_initial,f_noise_free,mu_opt,alpha_opt);
        norm(f_5-f_noise_free)/norm(f_noise_free)
        q = norm(B2'*f_5)%^2; %curl
        p = norm(B1*f_5)%^2; %div
        %(p -sqrt(p)*sqrt(q))/(p-q) - 1/(1+sqrt(q/p))
        par_opt(t) = p/(p+q)%(p -sqrt(p)*sqrt(q))/(p-q);  % which is equal to alpha = 1/(1+sqrt(q/p))
        alpha_opt = par_opt(t)%par_opt(t);
        %p^2/alpha_opt + q^2/(1-alpha_opt) == (p+q)^2
        %alpha = (p +sqrt(p)*sqrt(q))/(p-q); %the other quadratic solution
        %out of the bound
        if abs(alpha_opt-1) <1e-4
            alpha_opt = 1-1e-4;
        end
    end
    err_5_opt= err_5(end);
    if err_5_opt < best
        best = err_5_opt;
        mu_best_opt = mu_opt;
    end
    %cost_5 = 1/alpha_opt*norm(B1*f_5)^2 + 1/(1-alpha_opt)*norm(B2'*f_5)^2 +...
        %1/mu_opt*norm(f_5-f_initial)^2;
end
err_5_best = best;

close all
% M = [err_1];
% cHeader = {'err'};
% textHeader = strjoin(cHeader, ',');
% fid = fopen('example.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('example.csv',M,'-append','delimiter',',')