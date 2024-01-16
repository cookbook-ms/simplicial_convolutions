clear 
addpath('/Users/maosheng/OneDrive - Delft University of Technology/Topological filter/codes from 15 update/experiments/simplicial-regularization')

% this is a script to process flows on Minnesota transportation network
num_nodes = 2642; num_edges = 3303; num_tri = 53;
t1 = readtable('B1.csv');
t2 = readtable('B2t.csv');
t3 = readtable('coordinate.csv');
coordinate = t3{:,:}'; 
x_coord = coordinate(:,1);
y_coord = coordinate(:,2);
B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t';
num_edges= size(B1,2);
B1 = full(B1);
I = eye(num_edges);

L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;
% projectors;
P_g = B1'*pinv(B1*B1')*B1;
P_c = B2*pinv(B2'*B2)*B2';
P_h = eye(num_edges) - L1*pinv(L1);
%% spectrum
[u1l, lam_l]  = eig(L1l); eig_L1l = diag(lam_l); u_g = u1l(:,665:end); lam_g = eig_L1l(665:end); 
[u1, lam]= eig(L1); eig_L1 = diag(lam); u_h = u1(:,1:611);
[u1u, lam_u]= eig(L1u); eig_L1u = diag(lam_u); u_c = u1u(:,3251:end);lam_c = eig_L1u(3251:end);
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
f_noise_free = 5*f_harm + f_g + f_c ;
f_initial = 1*f_noise_free + 1*randn(num_edges,1);

flows_grad = P_g*f_initial; div = B1*f_initial; div_total = norm(div);
flows_harm = P_h*f_initial;
flows_curl = P_c*f_initial; curl = B2'*f_initial; curl_total = norm(curl);
err_initial = norm(f_initial-f_noise_free)/norm(f_initial);
%% study the spectrum of the underlying edge flow
f_tilde_noise_free = u1'*f_noise_free%f_noise_free;
figure;
plot(eig_L1,f_tilde_noise_free);
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
for alpha = 0.01:0.01:1
    if alpha==1
        continue;
    end
    for mu = 0.1:0.1:10
        f_1 = (I+ mu/alpha *L1l + mu/(1-alpha)*L1u)\f_initial;
        %f_1 = (I+ alpha *L1l + mu*L1u)\f_initial;
        err_1 = norm(f_1-f_noise_free)/norm(f_initial);
        if err_1 < best
            best = err_1;
            alpha_best = alpha;
            mu_best = mu;
        end
    end
end
err_1=best;
f_1 = (I+ mu_best/alpha_best *L1l + mu_best/(1-alpha_best)*L1u)\f_initial;

%%
best = 1e6;
for rho_2 = 0.01:0.01:10
    f_2 = (I+ rho_2 *L1l )\f_initial;
    err_2 = norm(f_2-f_noise_free)/norm(f_initial);
    if err_2 < best
        best = err_2;
        rho_2_best = rho_2;
    end
end
err_2=best;
f_2 = (I+ rho_2_best *L1l )\f_initial;
%%
best = 1e6;
for rho_3= 0.01:0.1:100
    f_3 = (I+ rho_3 *L1l + rho_3*L1u)\f_initial;
    err_3 = norm(f_3-f_noise_free)/norm(f_initial);
    if err_3 < best
        best = err_3;
        rho_3_best = rho_3;
    end
end
err_3=best;
f_3 = (I+ rho_3_best *L1l + rho_3_best*L1u)\f_initial;
%% MKL - optimal parameter design
best = 1e6;
%for mu = 0.1:0.1:10
mu = 10;
alpha = 0.4;
for t = 1:10 % number of iterations for alternating minimization
    f_5 = (I+ mu/alpha *L1l + mu/(1-alpha)*L1u)\f_initial;
    alpha = 1/(1+sqrt(norm(B2'*f_5)^2/norm(B1*f_5)^2));
end
err_5 = norm(f_5-f_noise_free)/norm(f_initial);
cost_5 = 1/alpha*norm(B2'*f_5)^2 + 1/(1-alpha)*norm(B1*f_5)^2;


%     if err_5 < best
%         best = err_5;
%         mu_2_best = mu;
%     end
%end
