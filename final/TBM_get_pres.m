function p_ratio = TBM_get_pres(M,beta,gamma)
%TBM_GET_PRES Summary of this function goes here
%   Detailed explanation goes here

p_ratio = 1 + 2.*gamma./(gamma+1) .* ( (M.^2).*(sin(beta).^2)-1 );

end

