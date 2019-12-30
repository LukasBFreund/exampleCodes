function x=cheby_scale_ab(z,a,b)
%%=========================================================================
% Scales Chebyshev nodes from [-1,+1] to [a,b]
% Inputs:  - n x 1 vector of Chebyshev nodes
%          - lower limit a
%          - upper limit b
% Output: n x 1 vector of scaled Chebyshev nodes
% See Elisa Faraglia Lecture Notes for PhD 21, Uni of Cambridge
% Lukas Freund, October 2018
%%=========================================================================

n = size(z,1); % dimension of nodes 
x = zeros(n,1); % storage object
x=((z+1)*(b-a))/2+a; %rescaling

end

