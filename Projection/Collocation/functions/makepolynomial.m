function mPoly = makepolynomial(vX1,vX2)
%%=========================================================================
% Creates polynomial combinations of inputs
%           - vX1: 1 x X1GridNum vector, row of Cheb. polynomials
%           - vX2: 1 x X2GridNum vector, row of Cheb. polynomials
% Output:   - mPoly: (max(X1GridNum,X2GridNum) *min(X1GridNum,X2GridNum))
% vector
% Lukas Freund, November 2018
%%=========================================================================

% Recover dimensions
X1GridNum = size(vX1,2);        % column numbers (order of Chebypol+1)
X2GridNum = size(vX2,2);        
XGridNumMax = max(X1GridNum,X2GridNum);
mPoly = [];

for iX2 = 1:X2GridNum
    for iX1 = 1:X1GridNum
        if (iX1+iX2-2)<=(XGridNumMax-1)
            mPoly = [mPoly kron(vX2(:,iX2), vX1(:,iX1))];
        end
end
end

end
