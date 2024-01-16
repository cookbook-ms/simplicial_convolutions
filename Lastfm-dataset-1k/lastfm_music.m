clear all
% this is a script to process flows on Chicago transportation network
addpath('~/Documents/chebfun')

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
 
%% creat a rating matrix -- scale from 1 to 5
% NOT SUCCESSFUL
% num_users = 10^4;
% rating_matrix = zeros(num_users, num_nodes);
% for i = 1:num_users
%     % set the number of ratings given by the i-th user
%     nnz_i = randi([20,30]);
%     for j = 1:nnz_i
%         % generate the random column indices, rated by the user i
%         nnz_j = randi(num_nodes,nnz_i,1);
%         rating_matrix(i,nnz_j(j)) = randi([1,5]);
%     end
% end
%
% spy(rating_matrix)
%
% % creat the edge flow matrix
% flow_mat = zeros(num_nodes);
% for i = 1:num_nodes
%     for j = 1:num_nodes
%         row_id = and(rating_matrix(:,i),rating_matrix(:,j));
%         if nnz(row_id) > 0
%             flow_mat(i,j) = sum(rating_matrix(row_id,i)-rating_matrix(row_id,j))/nnz(row_id);
%         end
%     end
% end
%% eigendecomposition for the SONG transition one
[U,Lam] = eig(L1); Lam = diag(Lam);
Lam(Lam(:)<1e-3) = 0; lam = uniquetol(Lam,0.003);
% the tolerance is the tol*abs(max)
% the maximal eigenvalue is around 10, we want the difference between
% 2.09-2.18 is zero, so we set tol = 0.1/36
% the harmonic space
U_H = U(:,1:1879);

[Ul,Lam_l] = eig(L1l); Lam_l = diag(Lam_l);
Lam_l(Lam_l(:)<1e-3) = 0; lam_l = uniquetol(Lam_l,0.1);
% the gradient space
U_G = Ul(:,2417:end);

[Uu,Lam_u] = eig(L1u); Lam_u = diag(Lam_u);
Lam_u(Lam_u(:)<1e-3) = 0; lam_u = uniquetol(Lam_u,0.01); % tol = 0.1/10
% the curl space
U_C = Uu(:,4971:end);

%% eigendecomposition for the ARTIST transition one
[U,Lam] = eig(L1); Lam = diag(Lam);
Lam(Lam(:)<1e-3) = 0; lam = uniquetol(Lam,0.006);
% the tolerance is the tol*abs(max)
% the maximal eigenvalue is around 10, we want the difference between
% 2.09-2.18 is zero, so we set tol = 0.1/64
% the harmonic space
U_H = U(:,1:477);

[Ul,Lam_l] = eig(L1l); Lam_l = diag(Lam_l);
Lam_l(Lam_l(:)<1e-3) = 0; lam_l = uniquetol(Lam_l,0.006);
% the gradient space
U_G = Ul(:,1342:end);

[Uu,Lam_u] = eig(L1u); Lam_u = diag(Lam_u);
Lam_u(Lam_u(:)<1e-3) = 0; lam_u = uniquetol(Lam_u,0.006);
% the curl space
U_C = Uu(:,1134:end);

%% flow spectrum
% analyze its frequency component
% these are the groundtruth subcomponent flows
f_h_tilde = U_H'*f; f_h = U_H*f_h_tilde;
f_g_tilde = U_G'*f; f_g = U_G*f_g_tilde;
f_c_tilde = U_C'*f; f_c = U_C*f_c_tilde;
f_tilde = U'*f;
% add some random divergence component in the groundtruth signal, the
% performance will be bad
% v = randi(10,num_nodes,1);
% f = f + B1'*v;
% add some known synthetic divergence component
% f_g_tilde_syn = zeros(num_edges,1);
% f_g_tilde_syn(1:5400) = 1;
% f_g_syn = U*f_g_tilde_syn; f = f+f_g_syn;
% check the divergence and curl
div = B1*f; curl = B2t*f;
t = pinv(B2'*B2)*curl;
t(abs(t)<5) = 0;
% observe the divergence component trend
figure(1);plot(f); figure(2);plot(curl);
figure(3);
stem(nonzeros(Lam_u),f_c_tilde);

% look for the large transition flows and the large curls, see if we can
% cluter the songs
% from curl with large or small values to the corresponding edge flows and
% three ndoes

% first study what the spectrum is like for the large curl signal leading
% to the edge flow
curl(abs(curl)<=5) = 0;
f_curl = B2*curl;
f_c_tilde_induced = U_C'*f_curl;
figure(4);
stem(nonzeros(Lam_u),f_c_tilde_induced);
hold on;
stem(nonzeros(Lam_u),f_c_tilde);
figure(5);
stem(t); hold on; stem(curl);
%% let us design a divergence-free, i.e., gradient removing filter, for the original dataset
L = 10;
Phi_G = [];
Phi_G = [Phi_G lam_l.^(0:1:L-1)];
% solve the LS problem
% for a gradient preserving filter, the parameter h = 1, preserving the rest
% alpha = pinv(0.2*eye(L)+Phi_G)*[1;zeros(length(lam_l)-1,1)];
% for the original data
g = [1;1;1;1;zeros(length(lam_l)-4,1)];
% this is for the synthetic case where we added some divergence component,
% we observe the gradient component trend
% g = [1;zeros(length(lam_l)-1,1)];
alpha = pinv(Phi_G)*g;
res1 = norm(Phi_G*alpha-g)/norm(g);
% build the filter
H1 = zeros(num_edges);
for l = 1:L
    H1 = H1 + alpha(l)*L1l^(l-1);
end

h_h1 = diag(U_H'*H1*U_H);
h_g1 = diag(U_G'*H1*U_G);
h_c1 = diag(U_C'*H1*U_C);

%% interpolation task
ratio = 0.001:0.1:0.999;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    for ii = 1:10
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
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
        f_unlabeled = pinv([B1*expanding_matrix;0.1*eye(nnz(mask_un))])...
            *[-B1*f_in;zeros(nnz(mask_un),1)];
        %f_unlabeled = lsqr([B1*expanding_matrix;0.1*eye(nnz(mask_un))],...
        %    [-B1*f_in;zeros(nnz(mask_un),1)]);
        % let us try the LSQR method
        f_jia2019 = f_in + expanding_matrix*f_unlabeled;
        corr_jia2019(i,ii) = corr(f,f_jia2019);
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out(i) = mean(corr_out_single(i,:),2);
    corr_jia2019_mean(i) = mean(corr_jia2019(i,:),2);
end

%% analyze its frequency component
f_h_tilde_o1 = U_H'*f_filtered1;
f_g_tilde_o1 = U_G'*f_filtered1;
f_c_tilde_o1 = U_C'*f_filtered1;
f_tilde_o1 = U'*f_filtered1;

% compare the divergence percentage change before and after the filter
div_in = norm(B1*f)/norm(f);
div_o = norm(B1*f_filtered1)/norm(f_filtered1);

%% new dataset
% consider a random flow on this music network, which are the transitions
% between each song or between each artist

%% PageRank to study the importance of some edges
% especially the edges with large upper degree
deg = diag(L1u);
[deg_max,pos] = maxk(deg,1); % pos tells us which edges have large degree
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% build the pagerank operator
k = 0.01;
H_P = inv(k*eye(num_edges) + L1_n);
% eigendecomposition of L1_n
[U_n,Lam_n,V_n] = eig(L1_n); % U_n is the right eigenvector; V_n' is the left one
Lam_n = diag(Lam_n);
Lam_n(Lam_n(:)<1e-3) = 0; lam_n = uniquetol(sort(Lam_n),0.1);
lam_n = linspace(0,1,21)';
% let us use the normalized Laplacian to formulate a filter to model the
% pagerank operator
% frequency response
h_p = 1./(k+lam_n);
for L = 1:11 % filter length
    Phi_n = lam_n.^(0:1:L-1);
    alpha = pinv(Phi_n) * h_p;
    % build the filter
    H_pr = zeros(num_edges);
    for l = 1:L
        H_pr = H_pr + alpha(l)*L1_n^(l-1);
    end
    err(L) = norm(H_pr - H_P)/norm(H_P); 
    % used for page ranking
    for i = 1:length(pos)
        ind_vec = zeros(num_edges,1);
        ind_vec(pos(i)) = 1;
        % the pagerank vector
        pg_vec(:,i) = H_P*ind_vec;
        % the filtering method
        pg_vec_H(:,i) = H_pr*ind_vec;
        figure;
        stem(pg_vec); hold on; stem(pg_vec_H);
        err_pr_vec(L,i) = norm(pg_vec - pg_vec_H)/norm(pg_vec);
    end
end


%% interpolation
% let us consider the normalzied Laplacian
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% the symmetrized version
L1_ns = D2^(-0.5)*L1_n*D2^(0.5);
ratio = 0.001:0.1:0.999;
for i = 1:length(ratio) 
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    %for ii = 1:1
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
        % the labeled
        f_in = f.*mask;
        % zero fill error
        corr_in(i) = norm(f-f_in)/norm(f);
        % build the sampling matrix
        sampling_mat = zeros(nnz(mask),num_edges);
        col_ind = find(mask==1);
        for j = 1:length(col_ind)
            sampling_mat(j,col_ind(j)) = 1;
        end
        
        % the filter method
        % the basic filter form
        % data-driven, design the filter coefficients
        F1 = [];
        for ll = 1:15
            F1 = [F1 L1^ll*f_in];
            % regularized LS to compute the filter coefficients
            h = inv(F1'*sampling_mat'*sampling_mat*F1 + 0.5*eye(ll))...
                *(F1'*sampling_mat'*sampling_mat)*f;
            % build the filter
            H = zeros(num_edges);
            for lll = 1:ll
                H = H + h(lll)*L1^lll;
            end
            f_1 = H*f_in;
            corr_filter_1(i,ll) = norm(f-f_1)/norm(f);
        end
end


%% heat diffusion operator approximation
t = 0.5;
H_hd = expm(-t*L1_n); 
h_hd = diag(U_n'*H_hd*U_n); 
figure; plot(1:length(Lam),h_hd);
% use the FIR filter to approximate this operator
L = 15; % the maximal order
H_hd_approx = zeros(size(H_hd));
for l = 0:L
    H_hd_approx  = H_hd_approx + (-t)^l/factorial(l)* L1_n^l;
    err_hd(l+1) = norm(H_hd-H_hd_approx)/norm(H_hd);
end

%% interpolation data driven
% let us consider the normalzied Laplacian
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% the symmetrized version
L1_ns = D2^(-0.5)*L1_n*D2^(0.5);
ratio = 0.1:0.1:0.8;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    %for ii = 1:1
    mask = zeros(num_edges,1);
    mask(randperm(numel(mask), M)) = 1;
    % the labeled
    f_in = f.*mask;
    % zero fill error
    corr_in(i) = norm(f-f_in)/norm(f);
    % build the sampling matrix
    sampling_mat = zeros(nnz(mask),num_edges);
    col_ind = find(mask==1);
    for j = 1:length(col_ind)
        sampling_mat(j,col_ind(j)) = 1;
    end
    sampling_mat = eye(num_edges);
    
    % the filter method
    % the basic filter form
    % data-driven, design the filter coefficients
    F1 = [];
    for ll = 1:10
        F1 = [F1 L1^(ll-1)*f_in];
        % regularized LS to compute the filter coefficients
        h = inv(F1'*sampling_mat'*sampling_mat*F1 + 0.1*eye(ll) ...
            + 0.5* F1'*B1'*B1*F1)...
            *(F1'*sampling_mat'*sampling_mat)*f;
        % build the filter
        H = zeros(num_edges);
        for lll = 1:ll
            H = H + h(lll)*L1^(lll-1);
        end
        f_1 = H*f_in;
        corr_filter_1(i,ll) = norm(f-f_1)/norm(f);
    end
    
%     % the subspace-varying filter
%     F2 = [];
%     for ll1 = 1:6
%         for ll2 = 1:5
%             F2 = [L1l^(0:1:ll1-1)*f_in L1u^(0:1:ll2)*f_in];
%             % regularized LS to compute the filter coeffieicents
%             h2 = inv(F2'*F2 + 0.1*eye(ll1+ll2))*F2'*f;
%             % build the filter
%             H2 = zeros(num_edges);
%             for lll1 = 1:ll1
%                 H2 = H2 + h2(lll1)*L1l^(lll1-1);
%                 
%             end
%             for lll2 = 1:ll2
%                 H2 = H2 + h2(ll1+lll2)*L1u^(lll2);
%             end
%             f_2 = H*f_in;
%             corr_filter_2(i,ll1,ll2) = norm(f-f_2)/norm(f);
%         end
%     end
end



