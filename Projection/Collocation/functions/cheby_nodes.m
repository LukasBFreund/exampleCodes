function z=cheby_nodes(n)
%%=========================================================================
% Function creates n Chebyshev nodes
% Inputs: n, the desired number of nodes 
% Output: n x 1 vector of Chebyshev nodes
% See http://www.wouterdenhaan.com/numerical/chebnode.m
% Lukas Freund, October 2018
%%=========================================================================

r	= max(1,n);
n	= (1:n)';
z	= cos((pi*(2*n-1))/(2*r));  %cos( (pi*(n-0.5))/r );

end
