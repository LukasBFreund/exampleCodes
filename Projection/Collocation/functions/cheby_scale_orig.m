function z=cheby_scale_orig(x,a,b)
%%=========================================================================
% Rescales Chebyshev nodes from [a,b] to [-1,+1]
% Inputs:  - n x 1 vector of scaled Chebyshev nodes
%          - lower limit a
%          - upper limit b
% Output: n x 1 vector of scaled Chebyshev nodes
% See Elisa Faraglia Lecture Notes for PhD 21, Uni of Cambridge
% Lukas Freund, October 2018
%%=========================================================================

n = size(x,1); % dimension of nodes 
z = zeros(n,1); % storage object
z = (2*(x-a))/(b-a)-1; 

end

