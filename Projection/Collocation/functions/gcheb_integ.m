function int_approx = gcheb_integ(f,x_l,x_u,m)
%%=========================================================================
% Gauss-Chebyshev quadrature (numerical integration) of a univariate function
% Inputs: 
%           - f: function, the integrand
%           - x_l: scalar, lower limit of integration
%           - x_u: scalar, upper limit of integration
%           - m: scalar, the desired order of integration 
% Output(s):
%           - int_approx: scalar, approximated integration
% Lukas Freund, October 2018
%%=========================================================================

%%Prepare nodes
cheb_nodes = cheby_nodes(m); % m x 1 vector of Chebyshev nodes (roots to Chebyshev polynomial of order m) defined on [-1,+1]
x_nodes = cheby_scale_ab(cheb_nodes,x_l,x_u); %scaled to [x_l,x_u]

%%Apply the approximation formula for Chebyshev quadrature
int_approx_arg1 = ((pi*(x_u-x_l))/(2*m));
int_approx_arg2 = sum(f(x_nodes).*((1-cheb_nodes.^2).^(1/2)));
int_approx = int_approx_arg1 * int_approx_arg2;

end

