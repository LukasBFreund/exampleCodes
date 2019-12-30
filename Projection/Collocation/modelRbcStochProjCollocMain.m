%%=========================================================================
% Stochastic RBC Model without labor, reversible investment
% AR(1) Case
% Projection Method: Collocation
% Lukas Freund, University of Cambridge
% November 2018
%%=========================================================================
%% --------------------------Current Notes------------------------------------

%% --------------------------Model------------------------------------
% V(k) = max_{C,K'} (C^(siggma-1)-1)/(siggma-1)) + betta*V(K',Z')
% with
% C + K' = Z*K^(alppha) + (1-deltta)*K
% log(Z') = rho_z * log(Z) + eps_z';

%% --------------------------Housekeeping------------------------------------
clear;
close all;
clc;
TimeStart = tic;
addpath('functions')

%% --------------------------Specification ---------------------------------
% Declare global variables (avoids functions requiring numerous inputs)
global KPolOrder lKMin lKMax mKChebPol vKGrid
global ZPolOrder lZMin lZMax mZChebPol vlZGrid 
global GHNodesNum ZNodesNum KNodesNum mGHWeights vGHGrid

% Structural parameters
betta = 0.95;                  % discount rate 
siggma=1.5;                             % inverse EIS
alppha = 0.33;                          % capital share 
deltta = 0.1;                         % capital depreciation rate  
sd_z = 0.02;                         % standard deviation, technology shock; changed this to 0.02 relative to problem set
rho_z = 0.8;                            % autoregressive coefficient, technology shock
mu_z = 0;                               % mean of technology shock
vParam = [alppha betta deltta siggma mu_z rho_z]'; % 6 x 1 vector

% Steady state (analytic)
Zss = 1;
Kss = (((1/betta)+deltta-1)/alppha)^(1/(alppha-1)); 
Yss = Zss*Kss^alppha;
Css = Yss - deltta*Kss;
Invss = deltta*Kss;

% Algorithm parameters 
KPolOrder = 4;                          % degree of polynomial, capital
KNodesNum = KPolOrder+1;                % # of nodes, capital
KCoeffNum = (KPolOrder+1)*(ZPolOrder+1); % # of coefficients, capital
ZPolOrder = 2;                          % degree of polynomial, technology shock
ZNodesNum = ZPolOrder+1;                % # of nodes, technology shock
GHNodesNum = 12;                        % # of nodes, Gauss-Hermite quadrature

%% --------------------------Grid ---------------------------------

% Capital grid and polynomials
KGridwidth = 0.5;
lKMin = log((1-KGridwidth)*Kss);                    % lower bound, in logs
lKMax = log((1+KGridwidth)*Kss);                    % upper bound, in logs
vKRoots = cheby_nodes(KNodesNum);                   % roots, in [-1,1]   
vKRoots = sort(vKRoots);                            % sort from low to high
vKGrid = exp(cheby_scale_ab(vKRoots,lKMin,lKMax));  % original scale
mKChebPol = cheby_pol(vKRoots,KPolOrder);           % Chebyshev polynomials for K

% Technology shock process: gaussian distribution => will use Gauss–Hermite quadrature
[vGHGrid,mGHWeights] = hernodes(GHNodesNum);
vGHGrid = vGHGrid*sqrt(2)*sd_z;           % GH nodes, with mGHWeights as associated weights

% Productivity grid and polynomials
lZMin = (mu_z-vGHGrid(1));                          % lower bound, in logs
lZMax = (mu_z+vGHGrid(1));                          % lower bound, in logs
vZRoots = cheby_nodes(ZNodesNum);                   % roots, in [-1,+1]
vZRoots = sort(vZRoots);                            % sort from low to high
vlZGrid = cheby_scale_ab(vZRoots,lZMin,lZMax);      % original scale (logs)
mZChebPol = cheby_pol(vZRoots,ZPolOrder);           % Chebyshev polynomials for Z

%% --------------------------Initial conditions ---------------------------------
% Initial guess
% Can be obtained from a linear approximation of the model or via homotopy
vTheta0 = [ 
    0.7212
    0.0493
    0.0017
    0.6527
   -0.0209
   -0.0005
    0.0138
    0.0016
    0.0000
   -0.0000
   -0.0000
    0.0000
];
vTheta0 = vTheta0(:);

%% --------------------------Main loop ---------------------------------
vTheta=fsolve('auxRbcStochProjCollocResiduals',vTheta0,[],vParam);

%% --------------------------Policy Functions ---------------------------------
CoeffNum = length(vTheta);
KGridNum = 100;
ZGridNum = 10;
ZGridInc = (lZMax-lZMin)/(ZGridNum-1);
vlZGridFine = [lZMin:ZGridInc:lZMax]';
vZGridFine = exp(vlZGridFine);
vlZTransf = cheby_scale_11(vlZGridFine,lZMin,lZMax);
mZChebPol = cheby_pol(vlZTransf,ZPolOrder);

% Capital
KGridInc = (lKMax-lKMin)/(KGridNum-1);
vlKGridFine = [lKMin:KGridInc:lKMax]';
vKGridFine = exp(vlKGridFine);                      % note this is a new (finer) grid
vlKTrans = cheby_scale_11(vlKGridFine,lKMin,lKMax);
mKChebPol = cheby_pol(vlKTrans,KPolOrder);

% Initialize
mC = [];            % KGridNum x ZNodesNum matrix
mOutput = [];
mKPrime= [];
mInv = [];

% Compute consumption, output, Kprime and investment
% Each is a  KGridNum x ZNodesNum matrix 

for iZ = 1:ZGridNum
    mPhi = makepolynomial(mZChebPol(iZ,:),mKChebPol); % approximated object
    mC = [mC exp(mPhi*vTheta)];
    mOutput = [mOutput exp(vlZGridFine(iZ))*vKGridFine.^(alppha)];
    mKPrime = [mKPrime mOutput(:,iZ)+(1-deltta)*vKGridFine-mC(:,iZ)]; 
    mInv = [mInv exp(vlZGridFine(iZ))*vKGridFine.^(alppha)-mC(:,iZ)];
end

%% --------------------------Plots ---------------------------------

%Plot policy functions
figure('Name','Policy Functions')
subplot(1,2,1)
f1 = mesh(vlZGridFine,vKGridFine,mKPrime);
view(-50,12);
xlabel('Productivity today, Z','FontSize',12,'fontname','times')
ylabel('Capital today,K','FontSize',12,'fontname','times')
zlabel('Capital tomorrow, K''','FontSize',12,'fontname','times')
title('Policy function: Capital');
subplot(1,2,2)
f2 = mesh(vlZGridFine,vKGridFine,mC);
view(-50,12);
xlabel('Productivity today, Z','FontSize',12,'fontname','times')
ylabel('Capital today, K','FontSize',12,'fontname','times')
zlabel('Consumption today','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
title('Policy function: Consumption');

% Plot Capital Policy function 
figure('Name','Policy Functions: Capital')
p1 = plot(vKGridFine,mKPrime(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vKGridFine,mKPrime(:,round(ZNodesNum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vKGridFine,mKPrime(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold on
p4 = plot(vKGridFine,vKGridFine,'--','color','k');
hold off
xlabel('Capital today,k','FontSize',12,'fontname','times')
ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2,p3,p4],'Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)','$45^{\circ}$ line');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy function(s)');

% Plot Consumption Policy function 
figure('Name','Policy Functions: Consumption')
p1 = plot(vKGridFine,mC(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vKGridFine,mC(:,round(ZNodesNum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vKGridFine,mC(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold off
xlabel('Capital today,K','FontSize',12,'fontname','times')
ylabel('Consumption today, C''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2,p3],'Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy function(s)');

% Plot investment policy function 
figure('Name','Policy Function: Investment')
p1 = plot(vKGridFine,mInv(:,1),'-','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vKGridFine,mInv(:,round(ZNodesNum/2+0.1)),'-','linewidth',1.5,'color',[1 0.5 0]);
hold on 
p3 = plot(vKGridFine,mInv(:,end),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
hold off
xlabel('Capital today, K','FontSize',12,'fontname','times')
ylabel('Investment, I''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2,p3],'Policy Rule (Worst Z)','Policy Rule (Median Z)','Policy Rule (Best Z)');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy function(s)');


%% --------------------------Done--------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']); 
