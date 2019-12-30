%%=========================================================================
% RBC Model without labor, without irreversibility constraint
% Endogenous Grid Point Method with polynomials
% No use of change of variables, so still finding a root
% Lukas Freund, University of Cambridge
% November 2018
%%=========================================================================

%% Housekeeping
%----------------------------------------------------------------------------
clear;
close all;
clc;
TimeStart = tic;
addpath('functions')
rng('default')

%% 1. Parameters 
%----------------------------------------------------------------------------

% Structural Model Parameters 
betta = (1.04)^(-1/4);                           % discount rate 
alppha = 0.36;                           % capital share 
deltta =0.025;                           % capital depreciation rate  
sigma_z = 0.02;                            % standard deviation, technology shock; changed this to 0.02 relative to problem set
rho_z = 0.98;                            % autoregressive coefficient, technology shock
gamma = 1.5;                              % inverse EIS

% Utility function
if gamma == 1
    fUtility = @(C) log(C);
else
    fUtility = @(C) ((C.^(1-gamma)-1)/(1-gamma));
end

% General algorithm parameters
SmoothingParam = 1;                     % Smoothing parameter for (fixed-point) iteration on value fn 
ConvTol = 1e-6;
PolOrder = 3;               

%% 2. Steady State
%----------------------------------------------------------------------------

% Steady state (analytic)
Zss = 1; 
Kss = (((1/betta)+deltta-1)/alppha)^(1/(alppha-1)); 
Yss = Zss.*Kss.^alppha;
Invss = deltta*Kss;
Css = Yss - Invss;

%% 3. State Variables Grid 
%----------------------------------------------------------------------------

% Unidimensional grid for capital 
KGridwidth = 0.5;
KMin = (1-KGridwidth)*Kss;                    % lower bound
KMax = (1+KGridwidth)*Kss;                    % upper bound
KGridNum = 10;                                % # of grid points for capital
vKGrid = linspace(KMin,KMax,KGridNum)';       % Unidimensional grid for capital, KGridNum x 1 

% Undimensional grid for productivity
ZGridwidth = 0.1;
ZMin = (1-ZGridwidth)*Zss;
ZMax = (1+ZGridwidth)*Zss;
ZGridNum = 10;
vZGrid = linspace(ZMin,ZMax,ZGridNum)';

% Two-dimensional tensor product grid: (KGridNum x ZGridNum) x 2
mGrid = [];
for iZ = 1:ZGridNum
    mGrid = [mGrid; [vKGrid ones(KGridNum,1)*vZGrid(iZ,1)]]; %vZGrid(iZ)
end
% aGrid = reshape(mGrid,[KGridNum,2,ZGridNum]) would be my usual notation 

mGrid_EGM = mGrid;                         % grid for EGM method
GridNum = KGridNum*ZGridNum;               % # of points in the 2D-grid
vKTGrid = mGrid(:,1);                      % capital grid points in 2D-grid
vZ = mGrid(:,2);                           % productivity grid points in 2D-grid

%% 4. Gauss-Hermite Quadrature
%----------------------------------------------------------------------------

QuadNum = 5;                             % # of integration nodes in GH quadrature
ShockNum = 1;                            % # of stochastic shocks
[x,w]=hernodes(QuadNum);
vGHWeights = w./sqrt(pi);
vGHNodes = x*sqrt(2)*sigma_z;           % GH nodes, with mGHWeights as associated weights
GHNodesNum = QuadNum;
mZ1 = vZ.^rho_z*exp(vGHNodes');         % future shocks in all grid points and all integration nodes, GridNum x GHNodesNum

%% 5. Initial Guess for Value Function on the Grid
%----------------------------------------------------------------------------

PolOrderInit = 2;                             % initial guess for value function is constructed using ordinary pol. of degree 2 
CoeffNum = 6;                                 % number of parameters in 2nd order polynomial 
mX0 = Polynomial_2d(mGrid,PolOrderInit);      % matrix of explanatory variables X0 on grid; columns: basis fn's of pol's of degree PolOrderInit

vVfCoeff = zeros(CoeffNum,1);                 % initial guess for coefficients of value function 

vVf = mX0*vVfCoeff;                           % initialize value fn. on the grid

% Initialize 
vC = vZ.*vKTGrid.^(alppha)*(Css/Yss);           % initialize by assuming a constant fraction css/yss of the period output goes to consumption
vKPrime = (1-deltta)*vKTGrid + vZ.*vKTGrid.^(alppha) - vC;
                                              % implied capital 

% Loop 
Metric = inf;                             

while Metric > ConvTol    
    for iG = 1:GHNodesNum
        mX1 = Polynomial_2d([vKPrime, mZ1(:,iG)],PolOrderInit);
        mEV(:,iG) = mX1*vVfCoeff;
    end
    
    vVfNew = fUtility(vC) + betta*mEV*vGHWeights;
    
    vVfCoeffNew = mX0\vVfNew;               % new vector of coefficients
    
    vVfCoeff = SmoothingParam*vVfCoeffNew+(1-SmoothingParam)*vVfCoeff;
    
    Metric = max(abs(1-vVfNew./vVf));       % difference between new and old value fn's 
    
    vVf = vVfNew;                           % update value function
end

%% 6. Main loop
%----------------------------------------------------------------------------

% Could run this for different polynomial orders
% for D = 2:4 
% PolOrder = D;

tic;
mX0 = Polynomial_2d(mGrid,PolOrder);             % ordinary polynomial of degree PolOrder
mX0Der = Polynomial_deriv_2d(mGrid,PolOrder);    % derivative
vVfCoeff = mX0\vVf;                             % initial guess for coefficients
vKCoeff = mX0\vKPrime;                              % initial guess for policy fn.
vKOld = inf(size(mGrid,1),1);                   % initialize capital choices for checking convergence

Metric = inf;

while Metric>ConvTol
    vKPrime = mGrid_EGM(:,1);                       % fix grid for next-period capital
    
    % Compute derivative of value function in the integration nodes  
    for iG = 1:GHNodesNum
        mXDer1_EGM = Polynomial_deriv_2d([vKPrime, mZ1(:,iG)],PolOrder);
        mVDer1_EGM(:,iG) = mXDer1_EGM*vVfCoeff;
    end
    
    % Compute expected derivative of next-period value fn
    vEVDer1 = mVDer1_EGM*vGHWeights;
    
    % Would solve for labor here if included elastically
    vC = (betta*vEVDer1).^(-1/gamma);       % Compute consumption using Euler eqn.
      
    % Next, compute current capital
    % for now, get capital using inefficient root-finding (should use Carrol's change
    % of variables!!
     fCapital=@(K,C,Z,KPrime)  (1-deltta)*K+Z.*K.^(alppha) - ( C + KPrime);

     for iK = 1:numel(vKTGrid)
     vKNew(iK) = fzero(@(K) fCapital(K,vC(iK),vZ(iK),vKPrime(iK)),vKTGrid(iK));
     end 
       
   % Recompute value function using the Bellman equation 
    for iG = 1:GHNodesNum
        mX1 = Polynomial_2d([vKPrime, mZ1(:,iG)],PolOrder);
        mEV(:,iG) = mX1*vVfCoeff;
    end
    
    vVfNew = fUtility(vC) + betta*mEV*vGHWeights;                                                                                                                   
    
    mGrid(:,1) = vKNew;                           % update grid points for current capital 
    
    mX0 = Polynomial_2d(mGrid,PolOrder);        % polynomial on current state variables
    
    vVfCoeff_New = mX0\vVfNew;                  % new coefficients for value fn, candidate
        
    vVfCoeff = SmoothingParam * vVfCoeff_New + (1-SmoothingParam) * vVfCoeff; % update coefficients
    
    % Check convergence of next period capital
    
    Metric = max(abs(1-vKNew./vKOld));
    
    vKOld = vKNew;
end

% Construct value function

MetricZ = 1e10;                                 % initially, convergence criterion not satisfied
mX0 = Polynomial_2d(mGrid,PolOrder);

while MetricZ>1e-10
    for iG = 1:GHNodesNum
        mX1 = Polynomial_2d([vKPrime, mZ1(:,iG)],PolOrder);
        mEV(:,iG) = mX1*vVfCoeff;
    end
    
    vVfNew = fUtility(vC) + betta*mEV*vGHWeights;
    
    vVfCoeffNew = mX0\vVfNew;               % new vector of coefficients
    
    vVfCoeff = SmoothingParam*vVfCoeffNew+(1-SmoothingParam)*vVfCoeff;
    
    MetricZ = max(abs(1-vVfNew./vVf));       % difference between new and old value fn's 
    
    vVf = vVfNew;                           % update value function
    
end
CPU(PolOrder) = toc;                      % Store running time
mCoeff(1:1+PolOrder+PolOrder*(PolOrder+1)/2,PolOrder) = vVfCoeff; % Store the solution coefficients for this particular order

% end 

%% 7. Plot policy function
%----------------------------------------------------------------------------
% Note that here it is K that varies rather than K1
mK0 = reshape(vKNew,KGridNum,ZGridNum);
mK0 = fliplr(mK0);                          % reading from left to right: worse to better productivity

% Plot Capital Policy function 
figure('Name','Policy Functions')
p1 = plot(mK0(:,1),vKGrid,'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(mK0(:,round(ZGridNum/2+0.1)),vKGrid,'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(mK0(:,end),vKGrid,'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold on
p4 = plot(vKGrid,vKGrid,'--','color','k');
hold off
xlabel('Capital today,k','FontSize',12,'fontname','times')
ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2,p3,p4],'Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)','$45^{\circ}$ line');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy function(s)');
 
%% --------------------------Done--------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']); 

