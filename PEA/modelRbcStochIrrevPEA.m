%%=========================================================================
% RBC Model without labor, irreversible investment
% PEA
% Lukas Freund, University of Cambridge
% November 2018
%%=========================================================================

%% --------------------------Housekeeping------------------------------------
clear;
close all;
clc;
TimeStart = tic;
addpath('functions')

%% --------------------------Algorithm description--------------------------------
% Model
% V(k) = max_{C,K'} (C^(siggma-1)-1)/(siggma-1)) + betta*V(K',Z')
% with
% C + K' = Z*K^(alppha) + (1-deltta)*K
% K'>=(1-deltta)*L <=> Inv>0
% log(Z') = rho_z * log(Z) + eps_z';
% We parameterize the conditional expectation with parameter vector Theta
% betta*E(t) [lambda(t+1)*(alppha*z(t+1)*k(t+1)^(alppha-1)+1-deltta)]
% Guess (if PolOrder ==2, otherwise exclude 2nd order terms
% Phi(K,Z; Theta) = exp(Theta_0+Theta_1 log(K)+Theta_2 log(Z)+Theta_3 log(K)^2 + Theta_4 log(Z)^2 + Theta_5 log(K)log(Z)
% We have from the equilibrium equations, with lambda being the Lagrange multiplier
% lambda(Theta) = Phi(K(Theta),Z(Theta);Theta) 
% C(t;Theta) = lambda(t;Theta)^(-1/siggma)
% K(t+1;Theta) = Z(t)*K(t;Theta)^alppha - C(t;Theta)+(1-deltta)*K(t;Theta)
% Given a series for K, Z, C, Lambda, can compute 
% phi(T+1;Theta) == lambda(t+1;Theta)*(alppha*Z(t+1;Theta)*K(t+1;Theta)^(alppha-1)+1-deltta)
% Regression: log(phi(t+1;Theta) = Theta_0+Theta_1 log(K(t,Theta))+Theta_2 log(K(t,Theta)...
% Recover Theta_hat, update Theta(i+1)=gamma*Theta_hat + (1-gamma)*Theta(i)
% To deal with irreversibility constraint, proceed as suggested in Marcet and Lorenzoni [1999]:
% Compute sequences under assumption that constraint not binding (mu=0)
% Check whether under this assumption Inv>0. If this is the case, set mu=0
% Otherwise set K(t+1)=(1-deltta)*K(t), compute C from the resource constraint
% and find mu(t) from the Euler equation, effectively mu(t) = C(t)^(-siggma) - lambda(t)

%% --------------------------Specification ---------------------------------
% Choose whether irreversibility constraint is included
ConstrIdx = 1;                          % Included if =1, excluded otherwise

% Parameterization of model deliberately chosen s.t. irreversibility
% constraint binds with a relatively high frequency
betta = 0.95;                           % discount rate 
alppha = 0.3;                           % capital share 
deltta = 0.09;                           % capital depreciation rate  
sigma_z = 0.08;                        % standard deviation, technology shock; changed this to 0.02 relative to problem set
rho_z = 0.8;                            % autoregressive coefficient, technology shock
siggma=1;                               % inverse EIS
mu_z = 0;                               % mean of technology shock

% PEA parameters
PolOrder = 2;                           % order of the polynomial used for approximating the cond. expectation
SmoothingAdaptive = 1;                   % If this is set to 1, then smoothing adjusts adaptively depending on magnitude of metric
SmoothingParam = 0.25;                           % smoothing parameter
NumData = 20000;                        % number of data points to compute the regression
NumDataInit = 500;                      % number of burn-in data points
NumDataTot = NumData + NumDataInit;
vT0 = NumDataInit+1:NumDataTot-1;       % for period t variables
vT1 = NumDataInit+2:NumDataTot;         % for period t+1 variables

% Collect parameters 
vParam = [mu_z alppha betta deltta rho_z sigma_z siggma NumData NumDataInit];

% Utility function
if siggma == 1
    fUtility = @(C) log(C);
else
    fUtility = @(C) ((C.^(1-siggma)-1)/(1-siggma));
end

% Production function
fProduction = @(k,z) z.*k.^alppha;

% Steady state (analytic)
Zss = 1;
Kss = (((1/betta)+deltta-1)/alppha)^(1/(alppha-1)); 
Yss = fProduction(Kss,Zss);
Css = Yss - deltta*Kss;
Invss = deltta*Kss;

% Algorithm specs
Tolerance = 1e-6;                       % convergence parameter
Metric = 1;                             % initialize convergence metric
Its = 0;                                 
%% --------------------------PEA --------------------------

% Generate a random series for the shock, Z
rng('default')
vEps = sigma_z*randn(NumDataTot,1);
vZ = zeros(NumDataTot,1);
vZ(1) = exp(mu_z+vEps(1));
for iT = 2:NumDataTot
    vZ(iT) = exp(rho_z*log(vZ(iT-1)))+(1-rho_z)*mu_z+vEps(iT); % note this defines Z (not log(Z))
end

% Make initial guess for parameter of approximating function
if PolOrder == 1
vTheta0 = zeros(3,1);                            % simplest guess (zeros has worked better than ones, in my experience)   
elseif PolOrder == 2
 %   vTheta0 = peaoginit(vEps,vParam);            
  %  vTheta0 = RbcStoch_pea_init(vEps,vParam);   % solves model relying on a log–linear approximation
%   vTheta0 = zeros(6,1);                        % typically doesn't converge 
vTheta0 = [
    0.3987
   -0.4208
   -0.6866
   -0.0508
   -0.2375
    0.2475
   ]; % guess from simpler version of model/with smaller shock variance (homotopy)

end

% Main loop
TimePeaIn=tic;

while Metric>Tolerance
    % Generate series
    vK = zeros(NumDataTot+1,1);    
    vlambda = zeros(NumDataTot,1);
    vmu =zeros(NumDataTot,1);                    % multiplier on the irreversibility constraint
    vC = zeros(NumDataTot,1);
    vInv = zeros(NumDataTot,1);
    mX = zeros(NumDataTot,length(vTheta0));      % matrix collecting series for regressors in polynomial fn.
    vK(1) = Kss;
    
    for iT = 1:NumDataTot
        if PolOrder == 2
        mX(iT,:) = [1 log(vK(iT)) log(vZ(iT)) log(vK(iT))*log(vK(iT)) log(vZ(iT))*log(vZ(iT)) log(vK(iT))*log(vZ(iT))]; %regressors
        elseif PolOrder == 1
        mX(iT,:) = [1 log(vK(iT)) log(vZ(iT))];
        end
        vlambda(iT) = exp(mX(iT,:)*vTheta0);        % recover lambda, assuming vTheta0
        vC(iT) = vlambda(iT)^(-1/siggma);           % recover consumption
        vInv(iT) = vZ(iT)*vK(iT)^alppha -vC(iT);
       
    if ConstrIdx ==1                                 % Irreversibility constraint included
        if vInv(iT)>0                                % in this case, constraint doesn't bind
            vK(iT+1) = (1-deltta)*vK(iT)+vInv(iT);
            vmu(iT) = 0;
        else                                          % constraint binds
            vK(iT+1) = (1-deltta)*vK(iT);
            vC(iT) =vZ(iT)*vK(iT)^alppha;       % total output, since Inv=0
            vmu(iT) = vC(iT)^(-siggma)-vlambda(iT);
             vInv(iT) = 0;
        end
    else                                               % No irreversibility constraint included in model
         vK(iT+1) = (1-deltta)*vK(iT)+vInv(iT);
    end
            
    end
    vphi = betta*vlambda(vT1).*(alppha*vZ(vT1).*vK(vT1).^(alppha-1)+1-deltta)... % recover dependent variable value
        -(1-deltta)*vmu(vT1); %note additional part due to irreversibility 
    vThetaCand = mX(vT0,:)\log(vphi);   % non-linear regression
    
      if SmoothingAdaptive==1     % Set smoothing parameter adaptively
    MetricCand = max(abs(vThetaCand-vTheta0));
    if MetricCand >0.1
        SmoothingParam = 0.1;
    elseif 0.1 <= MetricCand < 0.05
        SmoothingParam = 0.2;
    elseif 0.05 <= MetricCand <0.01
        SmoothingParam = 0.5;
    elseif MetricCand < 0.01
        SmoothingParam = 1;
    end
      end 
    
    vTheta = SmoothingParam*vThetaCand+(1-SmoothingParam)*vTheta0; % smoothing
    Metric = max(abs(vTheta-vTheta0));
    vTheta0 = vTheta; % update
    disp(vTheta(:)')
    disp(sprintf('Iteration: %d\tConv. crit.: %g',Its,Metric))
    Its=Its+1;
end

TimePea=toc(TimePeaIn);
disp(['PFI algorithm run time was ',num2str(TimePea),' seconds']); 

% Calculate probability that irreversibility constraint counted (0 if ConstrIdx ==0)
ConstrProb = sum(vmu(vT0)~=0)/length(vT0);
disp(['The irreversibility constraint was binding with a relative frequency of ',num2str(ConstrProb)]); 

%% --------------------------Policy Functions -------------------------------------

% Plot capital policy function 
figure('Name','Policy Function: Capital')
p1 = plot(vK(vT0),vK(vT1),'.','linewidth',1.5,'color',[0 0.45 0.75]);
hold on
p2 = plot(vK(vT0),vK(vT0),'--','color','k');
hold off
xlabel('Capital today, K','FontSize',12,'fontname','times')
ylabel('Capital tomorrow, K''','FontSize',12,'fontname','times')
legend1 = legend([p1,p2],'Policy Rule','$45^{\circ}$ line');set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy Function: Capital');

% Plot consumption policy function 
figure('Name','Policy Function: Consumption')
p1 = plot(vK(vT0),vC(vT0),'.','linewidth',1.5,'color',[0 0.45 0.75]);
hold off
xlabel('Capital today, K','FontSize',12,'fontname','times')
ylabel('Consumption today, C','FontSize',12,'fontname','times')
legend1 = legend([p1],'Policy Rule');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy Function: Consumption');

% Plot investment policy function 
figure('Name','Policy Function: Investment')
p1 = plot(vK(vT0),vInv(vT0),'.','linewidth',1.5,'color',[0 0.45 0.75]);
hold off
xlabel('Capital today, K','FontSize',12,'fontname','times')
ylabel('Investment today, Inv','FontSize',12,'fontname','times')
legend1 = legend([p1],'Policy Rule');set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex')
set(gca,'FontSize',12,'fontname','times')
title('Policy Function: Investment');

% Plot distribution of investment
figure('Name','Distribution of Investment')
hist(vInv, 100)
xlabel('Investment today, Inv','FontSize',12,'fontname','times')
ylabel('Frequency','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
title('Histogram: Investment');

% Plot typical path
TPlot=1000:1500-1;
figure('Name','Typical paths')
subplot(2,1,1)
p1 = plot(vInv(TPlot),'-','linewidth',1.5,'color','[0 0.45 0.75]');
hold on
p2 = plot(ones(length(TPlot),1)*Invss,'--','linewidth',1,'color','k');
legend1 = legend([p1 p2],'Investment','Steady state investment');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex');
xlabel('Time Period','FontSize',12,'fontname','times')
ylabel('Investment, Inv','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
title('Simulation: Investment');

subplot(2,1,2)
p1 = plot(vmu(TPlot),'-','linewidth',1.5,'color','[0 0.45 0.75]');
legend1 = legend([p1],'Multiplier, $\mu$');
set(legend1,'fontname','times','Location','best','FontSize',12,'interpreter','latex');
xlabel('Time Period','FontSize',12,'fontname','times')
ylabel('Multiplier, Inv','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
title('Simulation: Multiplier on Irreversibility Constraint');


%% --------------------------Done--------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']); 
