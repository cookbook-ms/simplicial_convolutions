clear all
% this is a script to process flows on Chicago transportation network
addpath('/Users/maosheng/Desktop/Topological filter/codes from 15 update/experiments')
t1 = readtable('B1.csv');
t2 = readtable('B2t.csv');
t3 = readtable('coordinate.csv');
t4 = readtable('flow_vec.csv');
num_nodes = 546; num_edges = 1088; num_tri = 112;

coordinate = t3{:,:}';
B1 = t1{:,:}; B2t = t2{:,:}; B2 = B2t';
L1l = B1'*B1; L1u = B2*B2t;
L1 = L1l + L1u;
f = t4{:,:};
f = f - mean(abs(f));

%%
p1 = figure;
plt_network(B1,f);
axis tight
exportgraphics(p1,'chicago_network_illustration.eps','BackgroundColor','none')



%% eigendecomposition
[U,Lam] = eig(L1); Lam = diag(Lam);
Lam(Lam(:)<1e-3) = 0; lam = uniquetol(Lam,0.06);
% the tolerance is the tol*abs(max)
% the maximal eigenvalue is around 10, we want the difference between
% 2.09-2.18 is zero, so we set tol = 0.1/10
% the harmonic space
U_H = U(:,1:431); lam_h = Lam(1:431);

[Ul,Lam_l] = eig(L1l); Lam_l = diag(Lam_l);
Lam_l(Lam_l(:)<1e-3) = 0; lam_l = uniquetol(Lam_l,0.06);
% the gradient space
U_G = Ul(:,544:end); lam_g = Lam_l(544:end);

[Uu,Lam_u] = eig(L1u); Lam_u = diag(Lam_u);
Lam_u(Lam_u(:)<1e-3) = 0; lam_u = uniquetol(Lam_u,0.06);
% the curl space
U_C = Uu(:,977:end); lam_c = Lam_u(977:end);
% note that in both the gradient frequency and the curl frequency. there
% are frequencies like 3, which occurs at both gradeint and curl
% frequencies

%% flow generation
% analyze its frequency component
f_h_tilde = U_H'*f; f_h = U_H*f_h_tilde;
f_g_tilde = U_G'*f; f_g = U_G*f_g_tilde;
f_c_tilde = U_C'*f; f_c = U_C*f_c_tilde;
f_tilde = U'*f;
% check the divergence and curl
div = B1*f; curl = B2t*f;
f - (f_h + f_g + f_c);
%% filter design -- gradient-preserving filter
% we shall reduce the divergence, thus design a frequency decaying
% first we want to test if the subspace-varying filter is able to deal with
% the situation where, e.g. 3 is a gradient and a curl frequency
% first remove the gradient frequency, so
% consider the perfect filter design, so the length should be larger than
% the number of distinct gradient frequencies
% but that will cause bad conditioning number, so we try L is small
L = 10;
Phi_G = [];
Phi_G = [Phi_G lam_l.^(0:1:L-1)];
% solve the LS problem
% for a gradient preserving filter, the parameter h = 1, preserving the rest
% alpha = pinv(0.2*eye(L)+Phi_G)*[1;zeros(length(lam_l)-1,1)];
alpha = pinv(Phi_G)*[0;ones(length(lam_l)-1,1)];
% build the filter
H1 = zeros(num_edges);
for l = 1:L
    H1 = H1 + alpha(l)*L1l^(l-1);
end

h_h1 = diag(U_H'*H1*U_H);
h_g1 = diag(U_G'*H1*U_G);
h_c1 = diag(U_C'*H1*U_C);

f_filtered1 = H1*f;
err_1 = norm(f_filtered1-f_g)/norm(f_g);
% analyze its frequency component
f_h_tilde_o1 = U_H'*f_filtered1;
f_g_tilde_o1 = U_G'*f_filtered1;
f_c_tilde_o1 = U_C'*f_filtered1;
f_tilde_o1 = U'*f_filtered1;
f_tilde_o1_rest = U'*(f-f_filtered1);

% compare the divergence percentage change before and after the filter
div_in = norm(B1*f)/norm(f);
div_o = norm(B1*f_filtered1)/norm(f_filtered1);

%%
p2 = figure;
plt_network(B1,f_g);
axis tight
exportgraphics(p1,'chicago_network_grad_subcomp_illustration.eps','BackgroundColor','none')


%% filter design -- curl-preserving filter
% we shall reduce the divergence, thus design a frequency decaying
% first we want to test if the subspace-varying filter is able to deal with
% the situation where, e.g. 3 is a gradient and a curl frequency
% first remove the gradient frequency, so
% consider the perfect filter design, so the length should be larger than
% the number of distinct gradient frequencies
% but that will cause bad conditioning number, so we try L is small
L = 10;
Phi_C = [];
Phi_C = [Phi_C lam_u.^(0:1:L-1)];
% solve the LS problem
% for a gradient preserving filter, the parameter h = 1, preserving the rest
% alpha = pinv(0.2*eye(L)+Phi_G)*[1;zeros(length(lam_l)-1,1)];
alpha = pinv(Phi_C)*[0;ones(length(lam_u)-1,1)];
% build the filter
H2 = zeros(num_edges);
for l = 1:L
    H2 = H2 + alpha(l)*L1u^(l-1);
end

h_h2 = diag(U_H'*H2*U_H);
h_g2 = diag(U_G'*H2*U_G);
h_c2 = diag(U_C'*H2*U_C);

f_filtered2 = H2*f ;
err_2 = norm(f_c-f_filtered2)/norm(f_c);
% analyze its frequency component
f_h_tilde_o2 = U_H'*f_filtered2;
f_g_tilde_o2 = U_G'*f_filtered2;
f_c_tilde_o2 = U_C'*f_filtered2;
f_tilde_o2 = U'*f_filtered2;

% compare the divergence percentage change before and after the filter
curl_in = norm(B2t*f)/norm(f);
curl_o = norm(B2t*f_filtered2)/norm(f_filtered2);

%% consider the universal design
lam_G_max = floor(max(lam_l));
% grid into 100 points
num_universal = 10;
lam_G_universal = [linspace(0,lam_G_max,num_universal)]';
L = 11;
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

%% figures
figure(1);
stem(f_h_tilde);  hold on; stem(f_h_tilde_o1); hold on;

figure(2); % check the gradient component in the frequency domain after filter H1
scatter(lam_g,f_g_tilde); hold on; scatter(lam_g,f_g_tilde_o1); hold on;

figure(3); % check the curl component in the frequency domain after filter H2
scatter(lam_c,f_c_tilde); hold on; scatter(lam_c,f_c_tilde_o2); hold on;


figure(4); % study the frequency response of H1
subplot(2,1,1);
figure; 
scatter(lam_g,h_g1,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
hold on;
scatter(lam_c,h_c1,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],...
    'MarkerEdgeColor',[0.9290 0.6940 0.1250]);
hold on;
scatter(0,mean(h_h1),'rs','filled');
ylim([0,1.2]);
xlim([0,11]);
set(gca,'fontsize',12);
legend('grad response h_{G}(\lambda)','curl response h_{C}(\lambda)',...
    'harmonic response h_{H}(0)');

subplot(2,1,2);
scatter(lam_g,h_g2,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
hold on;
scatter(lam_c,h_c2,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],...
    'MarkerEdgeColor',[0.9290 0.6940 0.1250]);
hold on;
scatter(0,mean(h_h2),'rs','filled');
ylim([0,1.2]);
xlim([0,11]);
xlabel('Frequency (\lambda)');
ylabel('Frequency response');
set(gca,'fontsize',12);
legend('grad response h_{G}(\lambda)','curl response h_{C}(\lambda)',...
    'harmonic response h_{H}(0)');


figure(5);
subplot(2,1,1);
scatter(Lam(1:700),f_tilde(1:700),'filled'); hold on;
scatter(Lam(1:700),f_tilde_o1_rest(1:700),'filled');
ylim([-0.5*10^(4),0.5*10^(4)])
xlim([0,3])
xticks(uniquetol([0;lam_c],0.03))
xtickformat('%.1f')
ax = gca;
ax.YAxis.Exponent = 3;
xlabel('Frequency (\lambda)');
ylabel('SFT');
set(gca,'fontsize',12);
legend('Input traffic flow', 'Output of $$ \mathbf{I}-\mathbf{H}_{G} $$',...
    'Interpreter','latex')

subplot(2,1,2);
scatter(Lam(700:end),f_tilde(700:end),'filled'); hold on;
scatter(Lam(700:end),f_tilde_o1_rest(700:end),'filled');
ylim([-0.5*10^(4),0.5*10^(4)])
xlim([2.9,10.5])
xticks(uniquetol([0;lam_c],0.05))
xtickformat('%.1f')
ax = gca;
ax.YAxis.Exponent = 3;
xlabel('Frequency (\lambda)');
ylabel('SFT');
set(gca,'fontsize',12);
legend('Input traffic flow', 'Output of $$ \mathbf{I}-\mathbf{H}_{G} $$',...
    'Interpreter','latex')

figure(6); % study the frequency response of H based on universal design
stem(lam_g,h_g); hold on;  stem(lam_c,h_c);

%% interpolation
% For classification, we consider the following labeled signal
% f = sign(randn(num_edges,1));
% f = randn(num_edges,1);
% let us consider the normalzied Laplacian
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% the symmetrized version
L1_ns = D2^(-0.5)*L1_n*D2^(0.5);
ratio = 0.1:0.1:0.6;
ratio = 0.5;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    %for ii = 1:1
    mask = zeros(num_edges,1);
    mask(randperm(numel(mask), M)) = 1;
    % the labeled
    f_in = f.*mask + mean(f).*randn(num_edges,1);
    %f_in = (eye(num_edges)+0.3*L1 + 0.6*L1^2)*f.*mask;
    % zero fill error
    err_in(i) = norm(f-f_in)/norm(f);
    corr_in(i) = corr(f,f_in);
    % build the sampling matrix
    sampling_mat = zeros(nnz(mask),num_edges);
    col_ind = find(mask==1);
    for j = 1:length(col_ind)
        sampling_mat(j,col_ind(j)) = 1;
    end
    %sampling_mat = eye(num_edges);
    
    % the filter method
    % the basic filter form
    % data-driven, design the filter coefficients
    F1 = [];
    for ll = 1:12
        F1 = [F1 L1^(ll-1)*f_in];
        % regularized LS to compute the filter coefficients
        %         h = inv(F1'*sampling_mat'*sampling_mat*F1 ...
        %             + 0.1*eye(ll) + 0.2*F1'*B2*B2'*F1)...
        %             *(F1'*sampling_mat'*sampling_mat)*f;
        h = inv(F1'*sampling_mat'*sampling_mat*F1 ...
            + 0.5*eye(ll))...
            *(F1'*sampling_mat'*sampling_mat)*f;
        % build the filter
        H = zeros(num_edges);
        for lll = 1:ll
            H = H + h(lll)*L1^(lll-1);
        end
        f_1 = H*f_in;
        %         % for classification
        %         f_1 = sign(H*f_in);
        err_filter_1(i,ll) = norm(f-f_1)/norm(f);
        corr_filter_1(i,ll) = corr(f,f_1);
    end
    
    % the subspace-varying filter
    F2 = [];
    for ll1 = 1:8
        F2 = [F2 L1l^(ll1-1)*f_in];
        for ll2 = 1:5
            F2 = [F2 L1u^(ll2)*f_in];
            % regularized LS to compute the filter coeffieicents
            % h2 = inv(F2'*F2 + 0.1*eye(ll1+ll2)+ 0.2*F2'*B2*B2'*F2)*F2'*f;
            h2 = inv(F2'*F2 + 0.5*eye(ll1+ll2))*F2'*f_in;
            % build the filter
            H2 = zeros(num_edges);
            for lll1 = 1:ll1
                H2 = H2 + h2(lll1)*L1l^(lll1-1);
                
            end
            for lll2 = 1:ll2
                H2 = H2 + h2(ll1+lll2)*L1u^(lll2);
            end
            f_2 = H2*f_in;
            %             % for classification
            %             f_2 = sign(H2*f_in);
            err_filter_2(i,ll1,ll2) = norm(f-f_2)/norm(f);
            corr_filter_2(i,ll1,ll2) = corr(f,f_2);
        end
        F2 = F2(:,1:ll1);
    end
end

%% interpolation based on the method in Jia2019, but combined with the subcomponent extraction filters
% the strategy is to first extract the subcomponent, then use the
% regularization methods in jia2019, as this could fit the underlying
% assumptions that the paper relies on, the divergence free and curl-free
ratio = 0.001:0.05:0.999;
for i = 1:length(ratio)
    M = floor(num_edges*ratio(i)); % the number of nonzero nodes
    for ii = 1
        mask = zeros(num_edges,1);
        mask(randperm(numel(mask), M)) = 1;
        % the labeled
        f_in = f.*mask;
        corr_in_single(i,ii) = corr(f,f_in); %norm(f-f_in)/norm(f);
        %         % the filter method
        %         f_filtered1 = H1*f_in;
        %         corr_out_single(i,ii) = corr(f,f_filtered1);
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
        corr_jia2019(i,ii) = corr(f,f_jia2019);  %norm(f-f_jia2019)/norm(f);
        % consider the filtering method then the regularization approach
        % first split the gradient subcomponents, which is curl free
        % regularization weights
        % THIS IS IMPRACTICAL, AS WE NEED TO ACCESS THE TRUE SIGNAL f TO
        % GENERATE THE SUBCOMPONENTS THEN DO MASK
        mu = 0.05;
        f_g = H1*f; f_in_g = f_g.*mask;
        %curl_in_g = norm(B2t*f_in_g)/norm(f_in);
        % this is based on penalizing the curl
        f_unlabeled_g = pinv([mu*eye(nnz(mask_un));B2t*expanding_matrix])...
            *[zeros(nnz(mask_un),1);-B2t*f_in_g];
        % then the rest components are divergence-free
        f_c = H2 * (f - f_g); f_in_c = f_c.*mask;
        %div_in_c = norm(B1*f_in_c)/norm(f_in);
        f_unlabeled_c = pinv([B1*expanding_matrix;mu*eye(nnz(mask_un))])...
            *[-B1*f_in_c;zeros(nnz(mask_un),1)];
        % then the harmonic component
        f_h = f - f_g-f_c; f_in_h = f_h.*mask;
        %         curl_in_h = norm(B2t*f_in_h)/norm(f_in);
        %         div_in_h = norm(B1*f_in_h)/norm(f_in);
        f_unlabeled_h = pinv([B1*expanding_matrix;mu*eye(nnz(mask_un));B2t*expanding_matrix])...
            * [-B1*f_in_h;zeros(nnz(mask_un),1);-B2t*f_in_h];
        
        f_out_g = f_in_g + expanding_matrix*(f_unlabeled_g);
        f_out_c = f_in_c + expanding_matrix*(f_unlabeled_c);
        f_out_h = f_in_h + expanding_matrix*(f_unlabeled_h);
        f_out = f_out_g + f_out_c + f_out_h;
        %         curl_out_g = norm(B2t*expanding_matrix*f_unlabeled_g)/norm(f_in);
        %         div_out_c = norm(B1*expanding_matrix*f_unlabeled_c)/norm(f_in);
        %         curl_out_h = norm(B2t*expanding_matrix*f_unlabeled_h)/norm(f_in);
        %         div_out_h = norm(B1*expanding_matrix*f_unlabeled_h)/norm(f_in);
        corr_out_g(i,ii) = corr(f_g,f_out_g); %norm(f-f_out)/norm(f); %
        corr_out_c(i,ii) = corr(f_c,f_out_c);
        corr_out_h(i,ii) = corr(f_h,f_out_h);
        corr_out(i,ii) = corr(f,f_out);
        
    end
    corr_in(i) = mean(corr_in_single(i,:),2);
    corr_out_g(i) = mean(corr_out_g(i,:),2);
    corr_out_c(i) = mean(corr_out_c(i,:),2);
    corr_out_h(i) = mean(corr_out_h(i,:),2);
    corr_out(i) = mean(corr_out(i,:),2);
    corr_jia2019_mean(i) = mean(corr_jia2019(i,:),2);
end

figure(6);
plot(ratio,corr_in);
hold on;
plot(ratio,corr_out);
hold on;
plot(ratio,corr_jia2019_mean);

%% PageRank to study the importance of some edges
[L1_n,D1,D2,D3] = compute_normalized_hodge_laplacian(B1,B2);
% build the pagerank operator
k = 0.01;
H_P = inv(k*eye(num_edges) + L1_n);
% eigendecomposition of L1_n
[U_n,Lam_n,V_n] = eig(L1_n); % U_n is the right eigenvector; V_n' is the left one
Lam_n = diag(Lam_n);
Lam_n(Lam_n(:)<1e-3) = 0;
% approximation design
lam_n = uniquetol(sort(Lam_n),0.1);
% universal design
lam_n = linspace(0,1,21)';
% let us use the normalized Laplacian to formulate a filter to model the
% pagerank operator
% frequency response
h_p = 1./(k+lam_n);
for L = 1:10 % filter lengthind_norm_h
    Phi_n = lam_n.^(0:1:L-1);
    alpha = pinv(Phi_n) * h_p;
    % build the filter
    H_pr = zeros(num_edges);
    for l = 1:L
        H_pr = H_pr + alpha(l)*L1_n^(l-1);
    end
    err(L) = norm(H_pr - H_P)/norm(H_P);
end
plot(lam_n,h_p)
stem(Lam_n, real(diag(U_n'*H_pr*U_n)))

%% pagerank 1 -- study several major edges
% especially the edges with large upper degree
% deg = diag(L1u);
% % pos tells us which edges have large degree，we consider the first three
% % largest
% [deg_max,pos] = maxk(deg,15);
% can we choose the ones with large number of neighbours to study
for i = 1:num_edges
    num_nbr_l(i) = nnz(L1l(:,i));
end
[num_nbr_l_max,pos_nbr_l] = mink(num_nbr_l,10);

for i = 1:num_edges
    num_nbr_u(i) = nnz(L1u(:,i));
end
[num_nbr_u_max,pos_nbr_u] = maxk(num_nbr_u,10);

pos = pos_nbr_l; 

%pos = pos_nbr_u;


% 
% p = figure;
% plt_network_page_rank(B1,pos);
% axis tight
% exportgraphics(p,'chicago_network_illustration_page_rank1.eps','BackgroundColor','none')
%  
% used for page ranking
for i = 1:length(pos)
    ind_vec = zeros(num_edges,1);
    ind_vec(pos(i)) = 1;
    % the pagerank vector
    pg_vec(:,i) = H_P*ind_vec;
    % the filtering method
    pg_vec_H(:,i) = H_pr*ind_vec;
    err_pr_vec(L,i) = norm(pg_vec - pg_vec_H)/norm(pg_vec);
    % look for the indices with a large pagerank value larger than 5
    [pg_vec_max5(:,i), ind_pg(:,i)] = maxk(pg_vec(:,i),10,...
        'ComparisonMethod','abs');
    [pg_vec_H_max5(:,i), ind_pg_H(:,i)] = maxk(pg_vec(:,i),10,...
        'ComparisonMethod','abs');    
end

% find the connected nodes of these most influenced edges
for i = 1:length(pos)
    for j = 1:10
        ind_pg_node1(j,i) = find(B1(:,ind_pg(j,i))==1);
        ind_pg_node2(j,i) = find(B1(:,ind_pg(j,i))==-1);
    end
end

%% make figures
% present the page rank results of edge 254 and 700 and 231, i.e., pos[1,3,6]
% p = figure;
% pg_254 = ind_pg(:,1);
% color = 'r';
% label = '254';
% plt_network_page_rank(B1,pg_254,color,label);
% 
% hold on;
% pg_700 = ind_pg(:,3);
% color = 'g';
% label = '700';
% plt_network_page_rank(B1,pg_700,color,label);
% 
% hold on; 
% pg_231 = ind_pg(:,6);
% color = 'y';
% label = '231';
% plt_network_page_rank(B1,pg_231,color,label);
% 
% hold on; 
% pg_231 = ind_pg;
% color = 'y';
% label = '231';
% plt_network_page_rank(B1,pg_231,color,label);

p = figure;
color = string(['#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F']);
for i = 1:5
    label = string(ind_pg(1,i));
    plt_network_page_rank(B1,ind_pg(:,i),color(i),label);
    hold on;
end

axis tight
exportgraphics(p,'chicago_network_illustration_page_rank3.eps',...
    'BackgroundColor','none')
 

% figure;
% % plot for the resuls of the ones with large lower degree
% subplot(2,2,1);
% for i = 1:3
%     scatter(ind_pg(:,i),abs(pg_vec_max5(:,i))./abs(max(pg_vec_max5(:,i))),...
%         'filled','LineWidth',2); hold on;
% end
% legend('254','388','700')
% set(gca,'fontsize',12)
% ylim([0,1])
% 
% subplot(2,2,2);
% for i = 4:6
%     scatter(ind_pg(:,i),abs(pg_vec_max5(:,i))./abs(max(pg_vec_max5(:,i))),...
%         'filled','LineWidth',2); hold on;
% end
% legend('790','229','231')
% set(gca,'fontsize',12)
% ylim([0,1])
% 
% % plot for the resuls of the ones with large upper degree
% subplot(2,2,3);
% for i = 1:3
%     scatter(ind_pg(:,i),abs(pg_vec_max5(:,i))./abs(max(pg_vec_max5(:,i))),...
%         'filled','LineWidth',2); hold on;
% end
% legend('180','181','183')
% set(gca,'fontsize',12)
% ylim([0,1])
% 
% subplot(2,2,4)
% for i = 7:9
%     scatter(ind_pg(:,i),abs(pg_vec_max5(:,i))./abs(max(pg_vec_max5(:,i))),...
%         'filled','LineWidth',2); hold on;
% end
% legend('249','255','257')
% set(gca,'fontsize',12)
% ylim([0,1])
% 
% 
% xlabel('Edge Index')
% ylabel('PageRank value')
% set(gca,'fontsize',12)

%% pagerank 2 -- study the influences on the gradient space of all edges
[num_nodes, num_edges] = size(B1);
num_tri = size(B2,2);
D2 = max(diag(abs(B2)*ones(num_tri,1)),eye(num_edges));
D1 = 2*diag(abs(B1)*D2*ones(num_edges,1));
D3 = 1/3*eye(num_tri);
L1_n = D2*B1'/D1*B1 + B2*D3*B2'/D2;

% compute the normalized hodge laplacian decomposition
% L1_n_s1 = D2^(-0.5)*L1_n*D2^(0.5); 
L1_n_s = D2^(0.5)*B1'/D1*B1*D2^(0.5) + D2^(-0.5)*B2*D3*B2'*D2^(-0.5);
% L1_n_s1 - L1_n_s2;
issymmetric(L1_n_s)

[U_n,Lam_n,V_n] = eig(L1_n_s); Lam_n = diag(Lam_n);
Lam_n(Lam_n(:)<1e-3) = 0; 
% the harmonic space
U_H_n = U_n(:,1:431);

[Ul_n,Lam_l_n] = eig(D2^(0.5)*B1'*(D2^(0.5)*B1')'); Lam_l_n = diag(Lam_l_n);
Lam_l_n(Lam_l_n(:)<1e-3) = 0; 
% the gradient space
U_G_n = Ul_n(:,544:end); 

[Uu_n,Lam_u_n] = eig(D2^(-0.5)*B2*(D2^(-0.5)*B2)'); Lam_u_n = diag(Lam_u_n);
Lam_u_n(Lam_u_n(:)<1e-3) = 0; 
% the curl space
U_C_n = Uu_n(:,977:end); 
%%
for i = 1:num_edges
    ind_vec = zeros(num_edges,1);
    ind_vec(i) = 1;
    % the pagerank vector
    pg_vec(:,i) = H_P*ind_vec;
    % the filtering method
    pg_vec_H(:,i) = H_pr*ind_vec;
    % evaluate in an aggregated fashion through 2-norm
    pg_norm(i) = norm(pg_vec(:,i));
    pg_H_norm(i) = norm(pg_vec_H(:,i));
    %err_pr_vec(i) = norm(pg_vec - pg_vec_H)/norm(pg_vec);
    % decompose and extract the gradient component, for its norm
    % we use the apprxoimately designed gradient preserving filter 
    pg_vec_g(:,i) = U_G_n'*pg_vec(:,i);
    pg_vec_H_g(:,i) = U_G_n'*pg_vec_H(:,i);
    % evaluate the curl norm
    pg_vec_c(:,i) = U_C_n'*pg_vec(:,i);
    pg_vec_H_c(:,i) = U_C_n'*pg_vec_H(:,i);
    % evaluate the harmonic norm
    pg_vec_h(:,i) = U_H_n'*pg_vec(:,i);
    pg_vec_H_h(:,i) = U_H_n'*pg_vec_H(:,i);
    
    % the absolute
    pg_norm_g(i) = norm(pg_vec_g(:,i));
    pg_H_norm_g(i) = norm(pg_vec_H_g(:,i));
    % the absolute
    pg_norm_c(i) = norm(pg_vec_c(:,i));
    pg_H_norm_c(i) = norm(pg_vec_H_c(:,i));
    % the absolute
    pg_norm_h(i) = norm(pg_vec_h(:,i));
    pg_H_norm_h(i) = norm(pg_vec_H_h(:,i));
end

% we should evaluate the normalized norm, i.e., wrt the total norm
% evaluate the gradient norm

for i = 1:num_edges
    % the relative
    pg_norm_g_rlt(i) = norm(pg_vec_g(:,i))/pg_norm(i);
    pg_H_norm_g_rlt(i) = norm(pg_vec_H_g(:,i))/pg_H_norm(i);
    % the relative
    pg_norm_c_rlt(i) = norm(pg_vec_c(:,i))/pg_norm(i);
    pg_H_norm_c_rlt(i) = norm(pg_vec_H_c(:,i))/pg_H_norm(i);
    % the relative
    pg_norm_h_rlt(i) = norm(pg_vec_h(:,i))/pg_norm(i);
    pg_H_norm_h_rlt(i) = norm(pg_vec_H_h(:,i))/pg_H_norm(i);
end

[~,ind_norm] = maxk(pg_norm,5);
[~,ind_H_norm] = maxk(pg_H_norm,5);
[~,ind_norm_g] = maxk(pg_norm_g,5);
[~,ind_H_norm_g] = maxk(pg_H_norm_g,5);
[~,ind_norm_c] = maxk(pg_norm_c,5);
[~,ind_H_norm_c] = maxk(pg_H_norm_c,5);
[~,ind_norm_h] = maxk(pg_norm_h,5);
[~,ind_H_norm_h] = maxk(pg_H_norm_h,5);

[~,ind_norm_g_rlt] = maxk(pg_norm_g_rlt,5);
[~,ind_H_norm_g_rlt] = maxk(pg_H_norm_g_rlt,5);
[~,ind_norm_c_rlt] = maxk(pg_norm_c_rlt,5);
[~,ind_H_norm_c_rlt] = maxk(pg_H_norm_c_rlt,5);
[~,ind_norm_h_rlt] = maxk(pg_norm_h_rlt,5);
[~,ind_H_norm_h_rlt] = maxk(pg_H_norm_h_rlt,5);

% let us investigate the structure of these edges, if they contribute to a
% triangle
for i = 1:length(ind_norm)
    ind_norm_node1(i) = find(B1(:,ind_norm(i))==1);
    ind_norm_node2(i) = find(B1(:,ind_norm(i))==-1);
    ind_norm_tri(i) = nnz(B2t(:,ind_norm(i))) ; 
end
for i = 1:length(ind_norm_g)
    ind_norm_g_node1(i) = find(B1(:,ind_norm_g(i))==1);
    ind_norm_g_node2(i) = find(B1(:,ind_norm_g(i))==-1);
    ind_norm_g_tri(i) = nnz(B2t(:,ind_norm_g(i))) ; 
end
for i = 1:length(ind_norm_c)
    ind_norm_c_node1(i) = find(B1(:,ind_norm_c(i))==1);
    ind_norm_c_node2(i) = find(B1(:,ind_norm_c(i))==-1);
    ind_norm_c_tri(i) = nnz(B2t(:,ind_norm_c(i))) ;
end
for i = 1:length(ind_norm_h)
    ind_norm_h_node1(i) = find(B1(:,ind_norm_h(i))==1);
    ind_norm_h_node2(i) = find(B1(:,ind_norm_h(i))==-1);
    ind_norm_h_tri(i) = nnz(B2t(:,ind_norm_h(i))) ;
end

%% illustration
p = figure;
color = string(['k';'g';'b';'r']);
ind_pg_norm = [ind_H_norm' ind_H_norm_g' ind_H_norm_c' ind_H_norm_h'];
ind_pg_norm_rlt = [ind_H_norm' ind_H_norm_g_rlt' ind_H_norm_c_rlt' ind_H_norm_h_rlt'];

for i = 1:4
    label = '';
    plt_network_page_rank(B1,ind_pg_norm(:,i),color(i),label);
    hold on;
end
%legend('Total','Grad','Curl','Harmonic')

axis tight
exportgraphics(p,'chicago_network_illustration_page_rank_norm_1.eps',...
    'BackgroundColor','none')

%%
figure; 
subplot(2,1,1);
% for subplot 1
%scatter(ind_norm,pg_norm(ind_norm),'sk'); hold on;
scatter(ind_H_norm,pg_H_norm(ind_H_norm),'sk','filled'); hold on;

%scatter(ind_norm_g,pg_norm_g(ind_norm_g),'sg'); hold on;
scatter(ind_H_norm_g,pg_H_norm_g(ind_H_norm_g),'sg','filled'); hold on;

%scatter(ind_norm_c,pg_norm_c(ind_norm_c),'sb'); hold on;
scatter(ind_H_norm_c,pg_H_norm_c(ind_H_norm_c),'sb','filled'); hold on

%scatter(ind_norm_h,pg_norm_h(ind_norm_h),'sr'); hold on;
scatter(ind_H_norm_h,pg_H_norm_h(ind_H_norm_h),'sr','filled'); 

% legend('2-norm by $\mathbf{H}_{PR}$', ...
%     '2-norm by $\mathbf{H}_1(\mathbf{L}_{1,n})$', ...
%     'Grad norm by $\mathbf{H}_{PR}$', ...
%     'Grad norm by $\mathbf{H}_1(\mathbf{L}_{1,n})$', ...
%     'Curl norm by $\mathbf{H}_{PR}$', ...
%     'Curl norm by $\mathbf{H}_1(\mathbf{L}_{1,n})$', ...
%     'Harmonic norm by $\mathbf{H}_{PR}$', ...
%     'Harmonic norm by $\mathbf{H}_1(\mathbf{L}_{1,n})$', ...
% 'Interpreter','latex')
legend('Norm of PageRank vec', 'Norm of grad comp', 'Norm of curl comp', ...
    'Norm of harmonic comp');
xlim([0,1100])
set(gca,'fontsize',12);
xlabel('Edge Index');

subplot(2,1,2);
% for subplot 2
%scatter(ind_norm,pg_norm(ind_norm)./pg_norm(ind_norm),'sk'); hold on;
scatter(ind_H_norm,pg_H_norm(ind_H_norm)./pg_H_norm(ind_H_norm),'sk','filled'); hold on;
scatter(ind_H_norm_g_rlt,pg_H_norm_g_rlt(ind_H_norm_g_rlt),'sg','filled'); hold on;
scatter(ind_H_norm_c_rlt,pg_H_norm_c_rlt(ind_H_norm_c_rlt),'sb','filled'); hold on
scatter(ind_H_norm_h_rlt,pg_H_norm_h_rlt(ind_H_norm_h_rlt),'sr','filled'); 
legend('Relative norm of PageRank vec', 'Relative norm of grad comp',...
    'Relative norm of curl comp', ...
    'Relative norm of harmonic comp');
xlim([0,1100])
set(gca,'fontsize',12);
xlabel('Edge Index');
%%
for i = 1:length(ind_pg)
    scatter(ind_pg(i),pg_vec(ind_pg(i)),'filled','MarkerFaceColor','k',...
        'MarkerEdgeColor','k');
    hold on;
    scatter(ind_pg(i),pg_vec_H(ind_pg_H(i)),'filled','MarkerFaceColor','b',...
        'MarkerEdgeColor','b');
end
hold on;
%scatter(180,0,'*',)
ylim([-6,11])
xlim([0,1088])
xticks(ind_pg);
%set(gca,'xaxisLocation','top')
xlabel('Edge Index');
set(gca,'fontsize',12);
legend('PageRank vec $$\mathbf{\pi}$$ by $$\mathbf{H}_{PR}$$', ...
    'PageRank vec $$\mathbf{\pi}$$ approx. by $$\mathbf{H}(\mathbf{L}_{1,n})$$',...
    'Interpreter','latex')

figure;
plot(1:L,err_pr_vec);
