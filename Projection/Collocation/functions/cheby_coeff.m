function theta=cheby_coeff(x,y,m)
%%=========================================================================
% Function calculates the coefficeints of a Chebyshev polynomial
% Inputs: 
%   - x: n x 1 vector of Chebyshev nodes defined on [-1,+1]
%   - y: n x 1 vector of function values defined on [a,b]
%   - m: order of Chebyshev polynomial
% Output: m+1 x 1 vector of Chebyshev coefficients
% See Elisa Faraglia Lecture Notes for PhD 21, Uni of Cambridge
% Lukas Freund, October 2018
% Note: Relative to lecture notes, the below def seems to omit a (-) sign
%%=========================================================================

%%Setup
n = size(x,1);  % Recover # of nodes
chebypol = zeros(n,m+1);    %Storage object for Chebyshev polynomials of order m (m+1 b/c first column contains ones)          
theta = zeros(m+1,1);             %storage object for m+1 coefficients

%%Get the polynomials first
chebypol(:,1) = chebyshevT(0,x);    %First column contains ones
for k=1:m
chebypol(:,k+1)=chebyshevT(k,x);    %Use of inbuilt function
end

%%Recover the coefficients, cf. lecture notes page 41
theta(1) = (1/n)*sum(y);          %theta_0
for j=2:size(theta)
    theta(j)=(2/n)*y'*chebypol(:,j);    %Note that y needs to be transposed for vector multiplication with the Chebyshev polynomial evaluated at the nodes x  
end

end
