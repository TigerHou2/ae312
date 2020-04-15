function theta = TBM_get_theta(M,beta,gamma)
%TBM_GET_THETA Summary of this function goes here
%   Detailed explanation goes here

theta = atan ( 2 .* cot(beta) .* ( ( (M.^2).*(sin(beta).^2)-1 ) ./ ...
                                 ( (M.^2).*(gamma+cos(2.*beta))+2 ) ) );

end

