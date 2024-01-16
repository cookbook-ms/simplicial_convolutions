clear all
% this is a script to process flows on Minnesota transportation network
num_nodes = 2642; num_edges = 3303; num_tri = 53;
t1 = readtable('B1.csv');
t2 = readtable('B2t.csv');
t3 = readtable('coordinate.csv');

coordinate = t3{:,:}'; 
x_coord = coordinate(:,1);
y_coord = coordinate(:,2);

B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t';
L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;

%% eigendecomposition
[U,Lam] = eig(L1); Lam = diag(Lam);
Lam(Lam(:)<1e-3) = 0; lam = uniquetol(Lam,0.0145);
% the tolerance is the tol*abs(max)
% the maximal eigenvalue is 6.8796, we want the difference between
% 2.09-2.18 is zero, so we set tol = 0.1/6.8796
% the harmonic space
U_H = U(:,1:611);

[Ul,Lam_l] = eig(L1l); Lam_l = diag(Lam_l);
Lam_l(Lam_l(:)<1e-3) = 0; lam_l = uniquetol(Lam_l,0.06);
% the gradient space
U_G = Ul(:,665:end);

[Uu,Lam_u] = eig(L1u); Lam_u = diag(Lam_u);
Lam_u(Lam_u(:)<1e-3) = 0; lam_u = uniquetol(Lam_u,0.0145);
% the curl space
U_C = Uu(:,3251:end);
% note that in both the gradient frequency and the curl frequency. there
% are frequencies like 3, which occurs at both gradeint and curl
% frequencies

%% flow generation
f = randn(num_edges,1);
% analyze its frequency component
f_h_tilde = U_H'*f;
f_g_tilde = U_G'*f;
f_c_tilde = U_C'*f;
f_tilde = U'*f;
%% filter design
% we shall reduce the divergence, thus design a frequency decaying 
% first we want to test if the subspace-varying filter is able to deal with
% the situation where, e.g. 3 is a gradient and a curl frequency
% first remove the gradient frequency, so 
% consider the perfect filter design, so the length should be larger than
% the number of distinct gradient frequencies
% but that will cause bad conditioning number, so we try L is small
L = 10;
Phi_G = [];
for l = 1:L
    Phi_G = [Phi_G lam_l.^(l-1)];
end
% solve the LS problem
% for a gradient removing filter, the parameter h = 1, preserving the rest
% alpha = pinv(0.2*eye(L)+Phi_G)*[1;zeros(length(lam_l)-1,1)];
alpha = pinv(Phi_G)*[1;zeros(length(lam_l)-1,1)];
% build the filter
H = zeros(num_edges);
for l = 1:L
    H = H + alpha(l)*L1l^(l-1);
end

h_h = diag(U_H'*H*U_H);
h_g = diag(U_G'*H*U_G);
h_c = diag(U_C'*H*U_C);

f_filtered = H*f;
% analyze its frequency component
f_h_tilde_o = U_H'*f_filtered;
f_g_tilde_o = U_G'*f_filtered;
f_c_tilde_o = U_C'*f_filtered;
f_tilde_o = U'*f_filtered;

% compare the divergence percentage change before and after the filter
div_in = norm(B1*f)/norm(f);
div_o = norm(B1*f_filtered)/norm(f_filtered);

%% consider the universal design
lam_G_max = ceil(max(lam_l));
% grid into 100 points
num_universal = 50;
lam_G_universal = [linspace(0,lam_G_max,num_universal)]';
L = 15;
Phi_G = lam_G_universal.^(0:1:L-1);
% the frequency response, if we set the gradient frequency response as 0,
% this will cause bad design, we set them as the inverse of the frequency
g = 0.02./(lam_G_universal(2:end)+0.1);
g = zeros(num_universal-1,1);
% compute the filter coefficient
alpha = pinv(Phi_G)*[1;g];

% build the filter
H = zeros(num_edges);
for l = 1:L
    H = H + alpha(l)*L1l^(l-1);
end
h_h = diag(U_H'*H*U_H);
h_g = diag(U_G'*H*U_G);
h_c = diag(U_C'*H*U_C);

f_o = H*f;
% analyze its frequency component
f_h_tilde_o = U_H'*f_o;
f_g_tilde_o = U_G'*f_o;
f_c_tilde_o = U_C'*f_o;

% compare the divergence percentage change before and after the filter
div_in = norm(B1*f)/norm(f);
div_o = norm(B1*f_o)/norm(f_o);


%%
% let us consider the regularization method
% suppress the gradient component
I = eye(num_edges);
H_rlg = inv(I+0.5*L1l);
h_h = diag(U_H'*H_rlg*U_H);
h_g = diag(U_G'*H_rlg*U_G);
h_c = diag(U_C'*H_rlg*U_C);
f_o = H_rlg*f;
% analyze its frequency component
f_h_tilde_o = U_H'*f_o;
f_g_tilde_o = U_G'*f_o;
f_c_tilde_o = U_C'*f_o;