function nu = PM(M,gamma)
%PM Prandtl-Meyer function

nu = sqrt((gamma+1)/(gamma-1)) ...
  .* atan(sqrt((gamma-1)/(gamma+1).*(M.^2-1))) ...
   - atan(sqrt(M.^2-1));

end

