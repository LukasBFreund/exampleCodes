function T = chebypol(x,n)
% Function for calculating Chebyshev polynominal up to order n
% See Heer and Maussner equation 11.54
% INPUTS
%   x:  mx1 vector of Chebyshev nodes defined on [-1,1]
%   n:  scalar, defining the order of interest
% OUTPUT:
%   T:  mx(n+1) vector

% (c) Simon Lloyd, September 2015

m = size(x,1);      % Number of grid points in x

T = ones(m,n+1);    % Define empty array for T

T(:,2) = x;
for i = 2:n
    T(:,i+1) = 2.*x.*T(:,i) - T(:,i-1);
end

