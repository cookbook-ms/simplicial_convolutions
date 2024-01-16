clear all;
close all;
% universal filter design - with a given continuous frequency response
% function of a variable, frequency, we want to design the filter
% coefficients.

% in this script, we will conduct a experiment to denoise the edge flow
% by designing the filter to behave like the regularization denoising
% solution

%% build the topology, i.e., B1 and B2 incidence matrix
b1 = [-1 -1 -1 0 0 0 0 0 0 0;...
    1 0 0 -1 0 0 0 0 0 0 ;...
    0 1 0 1 -1 -1 0 0 0 0;...
    0 0 1 0 1 0 -1 0 0 0 ;...
    0 0 0 0 0 0 1 -1 -1 0;...
    0 0 0 0 0 1 0 1 0 -1;...
    0 0 0 0 0 0 0 0 1 1];

b2 = zeros(10,3);
b2(1,1) = 1; b2(2,1) = -1; b2(4,1) = 1;
b2(2,2) = 1; b2(3,2) = -1; b2(5,2) = 1;
b2(8,3) = 1; b2(9,3) = -1; b2(10,3) = 1;

% build hodge laplacian matrix
Le = b1'*b1; L1u = b2*b2'; L1 = Le + L1u;

% eigendecomposition of laplacians
[u1l ~]  = eig(Le); eig_Le = eig(Le);
[u1 ~]= eig(L1); eig_L1 = eig(L1);
[u1u ~]= eig(b2*b2'); eig_L1u = eig(L1u);
[u1h ~] = eig([u1l(:,1:4), u1u(:,1:7)]*[u1l(:,1:4), u1u(:,1:7)]')

% generate an edge from the spectral domain which has all components with 1
% only at the zero frequencies, but zeros at the divergence and curl
% frequencies
f0 = u1 * [1;zeros(9,1)];
f0_tilde = u1' * f0; % perform fourier transform
norm_f0 = norm(f0)
% add some noise for robustness
noise = randn(10,1);
f = f0 + 1*noise;
f_tilde = u1' * f; % perform fourier transform
n_tilde = u1' * noise
err_ori = norm(f-f0)/norm(f0)

% the smallest that theoretically one can achieve
err_bound = norm(u1*([1;zeros(9,1)].*u1'*f) - f0)/norm(f0)

%% give the frequency response requirement
% build the frequency response requirement in continuous domain
% defined in the end of the script

%% perfrom the denoising by regularization solution
mu = [1e-2, 5e-2, 1e-1, 2.5e-1, 5e-1, 1, 2.5, 5, 10, 25, 50, 100]';
mu = 0.5
I = eye(size(L1));
A_lg = abs(b1'*b1-2*I);
L_lg = diag(A_lg) - A_lg;
for i = 1:length(mu)  % the larger the mu is, the better the denoising performance 
    H_rgl = inv(I+mu(i)*L1);
    H_rgl_freq = diag(u1'*H_rgl*u1);
    f_est_r = H_rgl*f;
    err_r(i) = norm(f_est_r-f0)/norm(f0);

    % line graph method
    H_lg = inv(I+mu(i)*L_lg);
    f_est_lg = H_lg*f;
    err_lg(i) = norm(f_est_lg-f0)/norm(f0);
end
%% make the plot of the desired filter frequency response
figure;
plot(eig_L1,H_rgl_freq,'k-.*','LineWidth',2); hold on;
grid on;
%% universal filter design
% first, we use fast algorithm to obtain the smallest and largest eigenvalues
% power iteration algorithm
v = ones(size(L1,1),1);
for k = 1:50
    v = L1*v;
    v = v/norm(v);
end

lambda_max = mean((L1*v)./v);
lambda_min = 0; % this is a doubt, does all hodge laplacian has a 0 eigenvalue?
% otherwise, the smallest eigenvalue can be found by inverse power
% iteration, but which takes a cubic order complexity again, not fast at
% all

% uniformly sample certain number of values from the eigenvalue interval
N = 10; % the number of samples

% eig_sampled = lambda_min + (lambda_max-lambda_min).*rand(N-2,1);
% % to make sure the smallest and largest eigenvalues are included in the
% % sampled eigenvalues
% eig_sampled = [lambda_min;eig_sampled;lambda_max];

% should we do uniform sampling or grid?
eig_sampled = [linspace(lambda_min,lambda_max,N)]';

% filter design
g = zeros(N,1);
A = [];
A_true = [];
L = 1:12; % filter lengths
h_coeff = [];
% to learn the regularization filter with topological filter 
for i = 1:N
    % generate the frequency response requirement vector
    g(i) = freq_response(eig_sampled(i));
end

g_true = 1./(1+mu.*eig_L1); % the true frequency response discrete 

% to remove the noise by topological filter to preserve the harmonic
% component, because in this situation we assume the ground truth flow have
% the small divergence and curl, Thus we design the topological filter
% universally from the simplicial domain and give the frequency response at
% 0 frequency to be 1 and 0 at the rest
% for i = 1:N
%     if eig_sampled(i) == 0 
%         g(i) = 1;
%     else
%         g(i) = 0;
%     end
% end


for l = L
    % build the system matrix
    A = [A, eig_sampled.^(l-1)];
    A_true = [A_true, eig_L1.^(l-1)];
    % solve the system equations by least squares solutions to obtain the
    % filter coefficients
    h_coeff = pinv(A)*g;
    h_coeff_true = pinv(A_true)*g_true;
    % build the topological filters
    H = zeros(size(L1));
    for i = 1:l
        H = H + h_coeff(i)*L1^(i-1);
    end
    
    H_true = zeros(size(L1));
    for i = 1:l 
        H_true = H_true + h_coeff_true(i)*L1^(i-1);
    end
    % perform the denoising
    f_est_tf = H * f;
    % compute the error
    err_tf(:,l) = norm(f_est_tf-f0)/norm(f0);
    % compare the frequency response designed with the given or the filter
    % comparison with the given filter
    err_filter(:,l) = norm(H - H_rgl); 
    err_filter_true(:,l) = norm(H-H_true);
    % the frequency response of the designed filter
    H_freq = diag(u1'*H*u1);
    if l==4 %l == 2 || l == 4 
        stem(eig_L1,H_freq,'-.d','lineWidth',2, 'MarkerSize',12); hold on;
    end
end
%%
xlabel('frequency ')
%ylabel('Filter frequency response')
set(gca,'fontsize',14)
legend('Desired filter','Grid design with L = 2',...
    'Grid design with L = 4')
xlim([0,6])
%% make some plots to show the errors
figure; % the Frobenius norm of the filters given and designed
plot(L,err_filter,'o-','LineWidth',2);
xlabel('Filter length (L)')
ylabel('$$ ||\mathbf{H}-\mathbf{G}||_F $$','Interpreter','latex')
title('Filter approximation by universal design')
grid on;
set(gca,'fontsize',12)
xlim([1 10])

figure;
plot(L,err_tf,'o-','LineWidth',2);
hold on;
plot(L,repmat(err_r,size(err_tf)),'*-','LineWidth',2);
legend('Topological filter','Regularization method')
xlabel('Filter length (L)')
ylabel('Normalized RMSE')
title('Denoising performance')
grid on;
set(gca,'fontsize',12)
xlim([1 10])

%% topological filter but designed in the spectral domain
A = []; % the system matrix to be updated
err_tf1 = [];
H_fr = [];
H_fr_leak = [];
a = [1;zeros(length(f0)-1,1)];
for l = L
    % build the system matrix
    A = [A, eig_L1.^(l-1)];
    h = pinv(A)*a;
    H = zeros(size(L1));
    % build the topological filters
    for i = 1:l
        H = H + h(i) * L1^(i-1);
    end
    % filtering process
    f_est = H * f;
    % compute the error % remember to change the exact component variable
    err_tf1(l) = norm(f_est - f0)/norm(f0);
    H_freq1 = diag(u1'*H*u1);
%     % check the filter frequency response
%     H_fr(:,L+1) = diag(u1' * H * u1);
%     % calculate the frequency leak
%     H_fr_leak(L+1) = norm(H_fr(:,L+1).* not(a));
end

%% plot the denoising performance
figure;
subplot(2,1,1);
plot(L,err_tf,'o-','LineWidth',2);
hold on;
plot(L,err_tf1,'o-','LineWidth',2);
xlabel('Filter length (L)')
ylabel('Normalized RMSE')
title('Topological filter method')
legend('Universal design','Spectral domain design')
grid on;
set(gca,'fontsize',12)
xlim([1 12])
ylim([0 1])

subplot(2,1,2);
plot(mu,err_r,'*-','LineWidth',2);
hold on;
plot(mu,err_lg,'*-','LineWidth',2);
title('Regularization method')
xlabel('Regularization weight \mu')
ylabel('Normalized RMSE')
legend('Hodge Laplacian rgl','Line-graph Laplacian rlg')
set(gca, 'XScale', 'log')
grid on;
set(gca,'fontsize',12)
xlim([min(mu) max(mu)])
ylim([0 3])
%% functions
function g = freq_response(x)
mu = 0.5
g = 1/(1+mu*x);
end
