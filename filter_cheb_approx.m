function H_cheb_approx_out = filter_cheb_approx(laplacian_operator,coeff,alpha,K_trnc)
coeff = coeff(1:K_trnc);
K = length(coeff);
for k = 1:K
    H_cheb_approx(k,:,:) = zeros(size(laplacian_operator));
    if k == 1
        H_cheb_approx(k,:,:) = eye(size(laplacian_operator));
    elseif k == 2
        H_cheb_approx(k,:,:) = 1/alpha*(laplacian_operator-alpha*eye(size(laplacian_operator)));
    else
        H_cheb_approx(k,:,:) = 2/alpha*(laplacian_operator-alpha*eye(size(laplacian_operator)))...
            *squeeze(H_cheb_approx(k-1,:,:)) - squeeze(H_cheb_approx(k-2,:,:));
    end
end
    H_cheb_approx_out = squeeze(sum(coeff.*H_cheb_approx,1));
end