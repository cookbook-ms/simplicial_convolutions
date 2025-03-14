function [estimate,error] = mkl_denoise(L1l,L1u,y,f0,mu,alpha)
% L1l, and L1u are the lower and upper laplacians
% y is the observation
% f0 is the true signal
% mu is the weight on the fitting term
% alpha is the weight controling the divergence and curl
I = eye(length(L1l));
estimate = (I+mu/alpha *L1l + mu/(1-alpha)*L1u)\y;
error = norm(estimate-f0)/norm(f0);

end