function [vRes] = auxRbcStochProjCollocResiduals(vTheta,vParam)
%%=========================================================================
% Computes vector of residuals (for collocation)
% Inputs:
%           - vTheta: vector of initial coefficients for the polynomial
%           - vParam: vector of structural parameters
% Output:   - vRes
% Lukas Freund, November 2018
%%=========================================================================
% Reference global variables
global KPolOrder lKMin lKMax vKGrid mKChebPol
global ZPolOrder lZMin lZMax vlZGrid mZChebPol
global GHNodesNum ZNodesNum KNodesNum mGHWeights vGHGrid

% Recover structurla parameters
alppha = vParam(1);
betta = vParam(2);
deltta = vParam(3);
siggma = vParam(4);
mu_z = vParam(5);
rho_z = vParam(6);

% Initialize
EulerRhs = [];                                          % right-hand side of Euler eqn.
EulerLhs = [];                                          % left-hand side of Euler eqn.
mX = [];                                                % Chebychev matrix 

%Check whether I can get this more directly via cheby fn. approximation!

for iK = 1:KNodesNum
    for iZ = 1:ZNodesNum
        vX0 = makepolynomial(mZChebPol(iZ,:),mKChebPol(iK,:)); 
        Ct = exp(vX0*vTheta);                                                     % C(t)
        mX = [mX;vX0];     
        % Compute K(t+1) 
        KPrime = exp(vlZGrid(iZ))*vKGrid(iK).^(alppha)+(1-deltta)*vKGrid(iK)-Ct;   % K(t+1)
        KPrimeTrans = cheby_scale_11(log(KPrime),lKMin,lKMax);                     % log(K(t+1)) in [-1;+1]
        vKPrimeCheb = cheby_pol(KPrimeTrans,KPolOrder);                            % Cheby. pol's for log(K_tr(t+1))
        % Compute Z(t+1)
        vZPrime = rho_z*vlZGrid(iZ)+(1-rho_z)*mu_z+vGHGrid;                        % Z(t+1), using GH-grid
        ZPrimeOrig = cheby_scale_11(vZPrime,lZMin,lZMax);                          % but that's not really scaled to [0,1], just same transf as K
        vZPrimeCheb = cheby_pol(ZPrimeOrig,ZPolOrder);                             % Cheby. pol's for Z(t+1)
        
        mX1 = makepolynomial(vZPrimeCheb,vKPrimeCheb);                             %Chebyshev polynomial matrix
        
        CPrime = exp(mX1*vTheta);                                                  % C(t+1)
              
        Integral = mGHWeights'*(betta*(alppha*exp(vZPrime)*KPrime.^(alppha-1)+1-deltta).*CPrime.^(-siggma))/sqrt(pi); %GH Integration
        EulerRhs = [EulerRhs;log(Integral)]; 
        EulerLhs = [EulerLhs;log(Ct.^(-siggma))];
         
    end
end

% Compute the residuals and transform into vector
vRes = (EulerLhs-EulerRhs);             % Residuals
vRes = vRes(:);

end
