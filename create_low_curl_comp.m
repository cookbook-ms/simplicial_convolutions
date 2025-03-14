function [f_c,f_c_tilde] = create_low_curl_comp(lam_c,mu_c,u_c)
% eig_L1l is the eigenvalues of the lower laplacian and ordered from small
% to big
% mu_g is the filtering parameter (I+mu*eig_L1l)^{-1}
% u_g is the gradient vectors
% I is the identity matrix
I = eye(length(lam_c));
f_c_tilde = (1+mu_c*lam_c).^(-1);
f_c = u_c*f_c_tilde;

end