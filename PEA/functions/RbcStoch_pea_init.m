function [vTheta] = RbcStoch_pea_NumDataInit(vEps,vParam)
%%=========================================================================
% Initializes parameters of stochastic growth model with a log-linear
% approximation. Note: no labor (but can be used for initialization)
% Inputs: 
%   - vEps: vector of shocks
%   - vParam: vector of model + algorithm parameters 
% Outputs:
%   - vTheta: vector of cofficients in approximating polynomial
% Implementation draws on code by F. Collard
% Lukas Freund, November 2018
%%=========================================================================

%% --------------------------Parameters: Algorithm-------------------------

% Parameters of algorithm 
NumCont = 3;                                % # of static equations (controls)
NumStateEndo = 1;                           % # of endogenous state v's
NumShock = 1;                               % # of shocks
NumS = NumStateEndo + NumShock;             % # of state variables
NumCostate = 1;                             % # of costate variables (aux/multiplier)
NumState = NumS+NumCostate;                 % # of state and costate variables

% Pre-allocate
mCC = zeros(NumCont,NumCont);
mCS = zeros(NumCont,NumState);
mSS0 = zeros(NumState,NumState);
mSS1 = zeros(NumState,NumState);
mSC0 = zeros(NumState,NumCont);
mSC1 = zeros(NumState,NumCont);
mSE = zeros(NumState,NumShock);

%% --------------------------Parameters: Economy---------------------------
% vParam = [mu_z alppha betta deltta rho_z sigma_z siggma NumData NumDataInit];
mu_z = vParam(1);
alppha = vParam(2);
betta = vParam(3);
deltta = vParam(4);
rho_z = vParam(5);
sigma_z = vParam(6);
siggma = vParam(7);

% And for algorithm
NumData = vParam(8);
NumDataInit = vParam(9);

%% --------------------------Steady State---------------------------------
Zss = 1;
Kss = (((1/betta)+deltta-1)/alppha)^(1/(alppha-1)); 
Yss =  Zss.*Kss.^alppha;
Invss = deltta*Kss;
Css = Yss - Invss;
Lambdass = Css^(-siggma);

%% --------------------------Matrix Coefficients---------------------------------
%% UNCLEAR
% Original has this, but it's not controls vs states
% y i 
% 1 2 
%
% k z c
% 1 2 3

% Output : Y = ZHat+ alppha*KHat
mCC(1,1) = 1;
mCS(1,1) = alppha;
mCS(1,2) = 1;

% Investment Inv = -Invss/Yss* ... -Css/Yss * ...
mCC(2,1) = 1;
mCC(2,2) = -Invss/Yss;
mCC(2,3) = -Css/Yss;

% Consumption: C-siggma
mCC(3,3) = -siggma;
mCS(3,3) = 1;

% Capital
mSS0(1,1) =1;
mSS1(1,1) = deltta-1;
mSC1(1,2) = deltta;

% Technology shock: -rho_z*...+eps
mSS0(2,2) = 1;
mSS1(2,2) = -rho_z;
mSE(2,1) = 1;

% Euler equation
mSS0(3,1) = (1-betta*(1-deltta));
mSS0(3,3) = -1;
mSS1(3,3) = 1;
mSC0(3,1) = (1-betta*(1-deltta));

%% --------------------------Solving the system---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                             %
%            Solving the system                               %
%                                                             %
%                                                             %
% X(t) = mCC^(-1)  mCS S(t)       %static                     %
%                                                             %
%                                                             %
% S(t+1)=(mSS0-mSC0 mCC^(-1)  mCS)(mSC1 mCC^(-1)mCS-mSS1)S(t) %
%                                                             %
%        +(mSS0-mSC0 mCC^(-1) mCS)E(t+1)                      %
%                                                             %
% S(t+1) = W S(t) + R E(t+1)                                  %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get W, which multipliers S(t)
M0=inv(mSS0-mSC0*inv(mCC)*mCS);
M1=(mSS1-mSC1*inv(mCC)*mCS);
W=-M0*M1;

[P,Mu] = eig(W);
MuAbs = abs(Mu);
[Mu,k] = sort(diag(MuAbs)); %sort diagonal entries from small to large
P=P(:,k);
Q=inv(P); % 3x3
 
% Direct solution
vGamma=-inv(Q(NumS+1:NumState,NumS+1:NumState))*Q(NumS+1:NumState,1:NumS);
mMSS=W(1:NumS,1:NumS)+W(1:NumS,NumS+1:NumState)*vGamma;
mPI=inv(mCC)*(mCS(:,1:NumS)+mCS(:,NumS+1:NumState)*vGamma); %for static equation
vMSE=[zeros(NumStateEndo,NumShock);eye(NumShock)]; %indicator where the shock is

S			= zeros(NumS,NumData+NumDataInit);
S(:,1)	    = vMSE*vEps(1);
for iT		= 2:NumData+NumDataInit
   S(:,iT)= mMSS*S(:,iT-1)+vMSE*vEps(iT);
end
% Generate series for the elements of the cond. expectation
lb			= vGamma*S;
lb			= Lambdass*exp(lb);
lb          = lb(:);
k			= log(Kss)+S(1,:);  %steady state plus shock component
k           = k(:);
ek			= exp(k);
z			= S(2,:);
z           = z(:);
ea			= exp(z);
T			= NumDataInit+1:NumDataInit+NumData-1;
T1			= NumDataInit+2:NumDataInit+NumData;

% Finally, run a regression
X			= [ones(NumData-1,1) k(T) z(T) k(T).*k(T) z(T).*z(T) k(T).*z(T)]; % regressors
y			= log(betta*lb(T1).*(alppha*ea(T1).*ek(T1).^(alppha-1)+1-deltta)); % regressand
vTheta		= X\y;

end
