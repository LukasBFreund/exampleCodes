function f_cheby_approx=cheby_fnapprox(f,z_l,z_u,n,m,N)
%%=========================================================================
% Wrapper for Chebyshev approximation of function f(z) defined on [z_l,z_u] using 
% Chebyshev polynomials of order m, n nodes, and a grid of N points
% Inputs: 
%   - f: function to be approximated
%   - z_l: lower limit of fn. domain
%   - z_u: upper limit of fn. domain
%   - n: number of Chebyshev nodes to be used, n
%   - m: order of Chebyshev polynomial
%   - N: number of grid points to be used for approximation (~finer grid)
% Output: 1 x N approximate function
% Lukas Freund, October 2018

%%=========================================================================

z_grid = linspace(z_l,z_u,N); %Create grid for approximation, not to be confused with # Chebyshev nodes 
f_z = f(z_grid);     % Grid evaluated using the original function 
x_grid = cheby_scale_orig(z_grid,z_l,z_u); %Scaled grid, to be used for function approximation  
z_nodes = cheby_nodes(n); % Find collocation points in [-1,+1] which are the roots of order m Chebyshev polynomials
x_nodes = cheby_scale_ab(z_nodes,z_l,z_u); %Transform these to corresponding grid points in [a, b]
y_nodes = f(x_nodes); %Use known information about f to get y_j=f(z_j)
theta = cheby_coeff(z_nodes,y_nodes,m); %Calculate the coefficients theta_i; note use of non-transformed nodes!
f_cheby_approx = cheby_approx(theta,x_grid,m);%Compute the approximation using the scaled grid   

end
