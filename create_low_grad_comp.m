function [f_g,f_g_tilde] = create_low_grad_comp(lam_g,mu_g,u_g)
% eig_L1l is the eigenvalues of the lower laplacian and ordered from small
% to big
% mu_g is the filtering parameter (I+mu*eig_L1l)^{-1}
% u_g is the gradient vectors
% I is the identity matrix
I = eye(length(lam_g));
f_g_tilde = (1+mu_g*lam_g).^(-1);
f_g = u_g*f_g_tilde;

end