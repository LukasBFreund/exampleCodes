%%=========================================================================
% RBC model without labor, reversible capital
% Envelope Condition Method (ECM) as per Maliar/Maliar (2013)
% Options for iterating on the value fn. or its derivative
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
betta = (1.04)^(-1/4);                  % discount rate 
alppha = 0.36;                          % capital share 
deltta = 0.025;                         % capital depreciation rate  
sigma_z = 0.02;                         % standard deviation, technology shock; changed this to 0.02 relative to problem set
rho_z = 0.98;                           % autoregressive coefficient, technology shock
gamma = 1.5;                            % inverse EIS
% A =  (1/betta-(1-deltta))/alppha;     % normalizing constant s.t. Kss = 1 
A = 1;                        

% Utility function
if gamma == 1
    fUtility = @(C) log(C);
else
    fUtility = @(C) ((C.^(1-gamma)-1)/(1-gamma));
end

% Production function
fProduction = @(K,Z) A*Z.*K.^(alppha);

% General algorithm parameters
ECMDVF = 1;                             % if ==0, iterature on value function, if ==1, iterate on derivative of VF
SmoothingParam = 1;                     % Smoothing parameter for (fixed-point) iteration on value fn 
ConvTol = 1e-6;
PolOrder = 4;               

%% 2. Steady State
%----------------------------------------------------------------------------

Zss = 1; 
Kss = ((((1/betta)+deltta-1)/alppha)/A)^(1/(alppha-1)); 
Yss = A*Zss.*Kss.^alppha;
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
ZGridwidth = 0.5;
ZMin = (1-ZGridwidth)*Zss;
ZMax = (1+ZGridwidth)*Zss;
ZGridNum = 10;
vZGrid = linspace(ZMin,ZMax,ZGridNum)';

% Two-dimensional tensor product grid: (KGridNum x ZGridNum) x 2
mTGrid = [];
for iZ = 1:ZGridNum
    mTGrid = [mTGrid; [vKGrid ones(KGridNum,1)*vZGrid(iZ,1)]]; %vZGrid(iZ)
end

TGridNum = KGridNum*ZGridNum;              % # of points in the 2D-grid
vKTGrid = mTGrid(:,1);                      % capital grid points in 2D-grid (T: tensor)
vZTGrid = mTGrid(:,2);                      % productivity grid points in 2D-grid (T: tensor)

%% 4. Gauss-Hermite Quadrature
%----------------------------------------------------------------------------

GHNodesNum = 5;
ShockNum = 1;                            % # of stochastic shocks
[x,w]=hernodes(GHNodesNum);
vGHWeights = w./sqrt(pi);
vGHNodes = x*sqrt(2)*sigma_z;           % GH nodes, with mGHWeights as associated weights
mZPrime = vZTGrid.^rho_z*exp(vGHNodes');     % future shocks in all grid points and all integration nodes, GridNum x GHNodesNum

%% 5. Initial Guess for Value Function on the Grid
%----------------------------------------------------------------------------

PolOrderInit = 2;                             % initial guess for value function is constructed using ordinary pol. of degree 2 
CoeffNum = 6;                                 % number of parameters in 2nd order polynomial 
mX = Polynomial_2d(mTGrid,PolOrderInit);      % matrix of explanatory variables X0 on grid; columns: basis fn's of pol's of degree PolOrderInit

vVfCoeff = zeros(CoeffNum,1);                 % initial guess for coefficients of value function 

vVf = mX*vVfCoeff;                           % initialize value fn. on the grid

% Initialize 
vC = A*vZTGrid.*vKTGrid.^(alppha)*(Css/Yss);  % initialize C by assuming a constant fraction Css/Yss
                                              % of period output goes to
                                              % consumption
vKPrime = (1-deltta)*vKTGrid + ...
           A*vZTGrid.*vKTGrid.^(alppha) - vC;  % implied K'

% Loop to initialize value function V' 
Metric = inf;                             

while Metric > ConvTol    
    for iG = 1:GHNodesNum
        mX1 = Polynomial_2d([vKPrime, mZPrime(:,iG)],PolOrderInit);
        mEV(:,iG) = mX1*vVfCoeff;
    end
    vVfNew = fUtility(vC) + betta*mEV*vGHWeights;
    vVfCoeffNew = mX\vVfNew;               % new vector of coefficients (A\B ~ inv(A)*B)
    vVfCoeff = SmoothingParam*vVfCoeffNew+(1-SmoothingParam)*vVfCoeff;
    Metric = max(abs(1-vVfNew./vVf));       % difference between new and old value fn's 
    vVf = vVfNew;                           % update value function
end

%% 6. Main loop: ECM iterating on value function
%----------------------------------------------------------------------------

% Could run this for different polynomial orders
% for D = 2:4 
% PolOrder = D;

mX = Polynomial_2d(mTGrid,PolOrder);                   % ordinary polynomial of degree PolOrder
mXDer = Polynomial_deriv_2d(mTGrid,PolOrder);          % derivative of polynomial
vVfCoeff = mX\vVf;                                     % initial guess for coefficients
vKCoeff = mX\vKPrime;                                  % initial guess for policy fn. (not needed for ECM)
vKPrimeOld = inf(size(mTGrid,1),1);                     % initialize K' for checking convergence
R = 1-deltta+A*alppha*vZTGrid.*vKTGrid.^(alppha-1);

Metric = inf;
ItsCount = 0;

TimeECMIn=tic;

while Metric > ConvTol
    vVDer = mXDer *vVfCoeff;                           % Derivative of value function, V_k (current)
    vMUtil = vVDer./(1-deltta+A*alppha*vZTGrid.*vKTGrid.^(alppha-1));
                                                        % U_c (marginal utility) from
                                                        % envelope condition (U_c=V_k/R)
    vC = vMUtil.^(-1/gamma);                               % Recover C
    vKPrime = fProduction(vKTGrid,vZTGrid)+(1-deltta)*vKTGrid-vC;
                                                        % Recover implied K'
                                                        
    if ECMDVF==0
        % Update value function                                      
        for iG = 1:GHNodesNum
            mX1 = Polynomial_2d([vKPrime, mZPrime(:,iG)],PolOrder);
            mEV(:,iG) = mX1*vVfCoeff;
        end
        vVfNew = fUtility(vC) + betta*mEV*vGHWeights;
        vVfCoeffNew = mX\vVfNew;               
    else    % Alternative: iterating on derivative of value function (ECM-DVF)
         for iG = 1:GHNodesNum
            mXPrimeDer = Polynomial_deriv_2d([vKPrime, mZPrime(:,iG)],PolOrder);
            vVDerPrime(:,iG) = mXPrimeDer*vVfCoeff;
         end
         vVDerNew = (1-deltta+A*alppha*vZTGrid.*vKTGrid.^(alppha-1)).*(betta*vVDerPrime*vGHWeights);
         warning('off')                         % system produces warning but LS problem still correctely processed
         vVfCoeffNew = mXDer\vVDerNew;
         
    end
    
     % Update parameters using smoothing
     vVfCoeff = SmoothingParam*vVfCoeffNew+(1-SmoothingParam)*vVfCoeff;

    % Check convergence of K'
    Metric = max(abs(1-vKPrime./vKPrimeOld));
    
    % Update
    vKPrimeOld = vKPrime;
    ItsCount = ItsCount + 1;
end

%% 7. Construct value function based on policy functions
%----------------------------------------------------------------------------

MetricZ = 1e10;                                 

while MetricZ>1e-10
    for iG = 1:GHNodesNum
        mX1 = Polynomial_2d([vKPrime, mZPrime(:,iG)],PolOrder);
        mEV(:,iG) = mX1*vVfCoeff;
    end
    
    vVfNew = fUtility(vC) + betta*mEV*vGHWeights;
    
    vVfCoeffNew = mX\vVfNew;               % new vector of coefficients
    
    vVfCoeff = SmoothingParam*vVfCoeffNew+(1-SmoothingParam)*vVfCoeff;
    
    MetricZ = max(abs(1-vVfNew./vVf));       % difference between new and old value fn's 
    
    vVf = vVfNew;                           % update value function
    
end
% CPU(PolOrder) = toc;                      % Store running time
% mCoeff(1:1+PolOrder+PolOrder*(PolOrder+1)/2,PolOrder) = vVfCoeff; % Store the solution coefficients for this particular order
TimeECM=toc(TimeECMIn);
disp(['ECM algorithm run time was ',num2str(TimeECM),' seconds']); 

%% 8. Illustrate results
%----------------------------------------------------------------------------

% Transform policy functions into matrix
mKPrime = reshape(vKPrime,KGridNum,ZGridNum);
mC = reshape(vC,KGridNum,ZGridNum);

%{
%Plot policy functions for capital and consumption in a mesh
figure('Name','Policy Functions')
subplot(1,2,1)
f1 = mesh(vZGrid,vKGrid,mKPrime);
view(-50,12);
xlabel('Productivity today, Z','FontSize',12,'fontname','times')
ylabel('Capital today,K','FontSize',12,'fontname','times')
zlabel('Capital tomorrow, K''','FontSize',12,'fontname','times')
title('Policy function: Capital');
subplot(1,2,2)
f2 = mesh(vZGrid,vKGrid,mC);
view(-50,12);
xlabel('Productivity today, Z','FontSize',12,'fontname','times')
ylabel('Capital today, K','FontSize',12,'fontname','times')
zlabel('Consumption today, C','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
title('Policy function: Consumption');

% Plot Capital Policy function 
figure('Name','Policy Functions')
p1 = plot(vKGrid,mKPrime(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vKGrid,mKPrime(:,round(ZGridNum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vKGrid,mKPrime(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold on
p4 = plot(vKGrid,vKGrid,'--','color','k');
hold off
xlabel('Capital today,k','FontSize',12,'fontname','times')
ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2,p3,p4],'Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)','$45^{\circ}$ line');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy function(s)');

%% 9. Simulation
%----------------------------------------------------------------------------

TSimPlot = 500-1;
TBurn = 200; % # of burn-in periods to eliminate the effect of initial conditions
TSim = 10000;
T = TSim+TBurn;

% Construct exogenous productivity series
vEpsSt = randn(T,1);                   % Generate standard normal values
vEps = sigma_z*vEpsSt;                 % standard deviation is sigma_z

vZSim(1,1) = Zss;                            % initial productivity is 1
for iT = 2:T
vZSim(iT,1) = vZSim(iT-1).^rho_z.*exp(vEps(iT));
end

% Simulate endogenous variables
vKSim(1,1) = Kss;                            % initial capital is at SS

for iT = 1:T
 VDer = Polynomial_deriv_2d([vKSim(iT,1),vZSim(iT,1)],PolOrder)*vVfCoeff; % V_k
 MUtil = VDer./(1-deltta+A*alppha*vZSim(iT,1).*vKSim(iT,1).^(alppha-1));  % MU_c
 vCSim(iT,1) = MUtil.^(-1/gamma);
 vYSim(iT,1) = A*vZSim(iT,1).*vKSim(iT,1).^(alppha);
 vInvSim(iT,1) = vYSim(iT,1) - vCSim(iT,1);
 vKSim(iT+1,1) = (1-deltta)*vKSim(iT,1)+vInvSim(iT,1);       
end

% Discard burn in values
vZSim = vZSim(TBurn+1:T,1);
vKSim = vKSim(TBurn+1:T,1);
vCSim = vCSim(TBurn+1:T,1);
vInvSim = vInvSim(TBurn+1:T);
vYSim = vYSim(TBurn+1:T);

% Compute moments
vSimMeans = mean([vKSim vInvSim vYSim vCSim])';
SimSD = sqrt(var([vKSim vInvSim vYSim vCSim]))';
S.Name={'Capital, k'; 'Investment, i'; 'Output, y'; 'Consumption, c'};
S.Mean=vSimMeans;
S.SD =SimSD;
disp('Moments of k, i, y and c, respectively:');
struct2table(S)

% Plot result (for last TSimPlot periods)
figure('Name','Simulation');
subplot(2,2,1);
plot(vZSim(end-TSimPlot:end),'linewidth',1.6,'color','k');
ylabel('Productivity, z','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,2);
plot(vYSim(end-TSimPlot:end),'linewidth',1.6,'color','k');
ylabel('Output, y','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,3);
plot(vInvSim(end-TSimPlot:end),'linewidth',1.6,'color','k');
ylabel('Investment, i','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,4);
plot(vCSim(end-TSimPlot:end),'linewidth',1.6,'color','k');
ylabel('Consumption, c','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
%}

%% --------------------------Done--------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']); 
