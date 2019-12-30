function f_approx=cheby_approx(theta,x,m)
%%=========================================================================
% Inputs: Approximates function f 
%   - m: order of Chebyshev polynomial
%   - x: N x 1 vector of (scaled) Chebyshev nodes, note that N is the number of grid points, not the number of Chebyshev nodes
%   - theta: (m+1) x 1 vector of Chebyshev coefficients
% Output: 1 x N vector approximation of functin f 
% Lukas Freund, October 2018
%%=========================================================================

%%Setup
N = size(x,2);  % Recover # of nodes
chebypol = zeros(N,m+1); 

%%Create the matrix of
chebypol(:,1) = chebyshevT(0,x);    %First column contains ones
for k=1:m
chebypol(:,k+1)=chebyshevT(k,x);    %Use of inbuilt function
end

%%Approximation
f_approx = theta'*chebypol'; %note that theta is (m+1) x 1 while chebypol is N x (m+1), so theta' * chebypol' is (1x(m+1))*(m+1)xN => 1x N
end
