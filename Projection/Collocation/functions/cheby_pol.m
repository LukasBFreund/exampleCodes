function vChebypol = cheby_pol(vChebynodes,ChebyOrder)
%%=========================================================================
% Calculate Chebyshev polynominal up to order ChebyOrder
% Inputs:
%           - vChebynodes: n x 1 vector, Chebyshev nodes defined on [-1,1]
%           - m: scalar, defining the order of interest
% Output:   - vChebypol, n x (ChebyOrder+1) vector of Chebyshev polynmials
% Lukas Freund, November 2018
%%=========================================================================

GridNum = size(vChebynodes,1);
vChebypol = ones(GridNum,ChebyOrder+1);

vChebypol(:,2) = vChebynodes;
for iN = 2:ChebyOrder
    vChebypol(:,iN+1) = 2.*vChebynodes.*vChebypol(:,iN) - vChebypol(:,iN-1);
end

end
