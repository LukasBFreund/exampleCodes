%%=========================================================================
% Solves (non-linear) baseline NK model of J. Gali (2015), ch. 3
% Using P. Rendahl's Linear Time Iteration (LTI) method (Rendahl, 2017)
% Implementation mostly follows J. Pfeifer's Dynare code for comparison
% However, interest rates and inflation are not annualized
% Lukas Freund, University of Cambridge
% Written with Matlab R2018a
% October 2018
%%=========================================================================

%% --------------------------Housekeeping------------------------------------
clear;
close all;
clc;
TimeStart = tic;
%% --------------------------Specification------------------------------------
Shock = 'MP'; %choose between TFP, MP, Preference shock
IrfHor = 11; %choose length of IRF
Tol = 1e-13; %choose tolerance level for LTI algorithm, default: 1e-13

%% --------------------------Parameterization------------------------------------
siggma = 1; %inverse EIS
varphi= 5 ; %inverse Frisch elasticity
phi_pi = 1.5; %inflation feedback Taylor rule
phi_y  = 0.125; %output feedback Taylor rule
theta = 3/4; %Calvo parameter
rho_z  = 0.5; %autocorrelation, preference shock
rho_a  = 0.9; %autocorrelation, technology shock
rho_nu = 0.5; %autocorrelation, interest rate shock
betta  = 0.99; %discount factor
alppha = 1/4; %capital share
epsilon = 9; %demand elasticity 
tau = 0;  %labor subsididy
sig_z = 0.005; %size of preference shock (% deviation)
sig_a = 0.01; %size of technology shock
sig_nu = 0.0025; %size of monetary policy shock

%% --------------------------Deterministic Steady State--------------------------------
% Define non-stochastic SS following Gali/Pfeiffer
A_ss=1; %productivity
Z_ss=1; %preference shifter
S_ss=1; %price dispersion
nu_ss=1; %monetary policy shock
Pi_star_ss=1; %optimal inflation level
P_ss=1; %price level
MC_ss=(epsilon-1)/epsilon/(1-tau); %marginal cost
R_ss=1/betta; %nominal interest rate 
realinterest_ss=R_ss; %real interest rate 
Pi_ss=1; %inflation level
Q_ss=1/R_ss; %bond price
N_ss=((1-alppha)*MC_ss)^(1/((1-siggma)*alppha+varphi+siggma)); %hours worked
C_ss=A_ss*N_ss^(1-alppha); %consumption
W_real_ss=C_ss^siggma*N_ss^varphi; %real wage 
Y_ss=C_ss; %output
x_aux_1_ss=C_ss^(-siggma)*Y_ss*MC_ss/(1-betta*theta*Pi_ss^(epsilon/(1-alppha))); %auxiliary price setting recursion 1
x_aux_2_ss=C_ss^(-siggma)*Y_ss/(1-betta*theta*Pi_ss^(epsilon-1)); %auxiliary price setting recursion 2

% Summarize SS and get number of variables
xss = [A_ss,Z_ss,S_ss,Pi_star_ss,P_ss,MC_ss,R_ss,realinterest_ss,Pi_ss,Q_ss,N_ss,C_ss,W_real_ss,Y_ss,x_aux_1_ss,...
    x_aux_2_ss,nu_ss];
NumVar = size(xss,2);

%% --------------------------Specify Stochastic System--------------------------------
% Specify symbolic variables 
syms A Z S Pi_star P MC R realinterest Pi Q N C W_real Y x_aux_1 x_aux_2 nu realinterest_ann
syms lA lZ lS lPi_star lP lMC lR lrealinterest lPi lQ lN lC lW_real lY lx_aux_1 lx_aux_2 lnu lrealinterest_ann
syms fA fZ fS fPi_star fP fMC fR frealinterest fPi fQ fN fC fW_real fY fx_aux_1 fx_aux_2 fnu frealinterest_ann

% Select the variables that will be shown in the plot and name them
VarPlot = [7,8,9,11,13,14]; 
VarNames = {'Nominal Interest Rate','Real Interest Rate','Inflation','Hours','Real Wage','Output'};

% Specify equations (symbolic) 
e1=W_real-(C.^siggma.*N.^varphi);  %FOC labor
e2=Q-(betta*(fC./C).^(-siggma).*(fZ./Z)./fPi); %Euler equation
e3=R-(1./Q); %definition nominal interest rate
e4=Y-(A.*(N./S).^(1-alppha)); %aggregate output
e5=R-(realinterest.*fPi); %definition real interest rate
e6=R-(1/betta.*Pi.^phi_pi.*(Y./Y_ss).^phi_y.*exp(nu)); %monetary policy rule
e7=C-Y; %market clearing
e8=-log(A)+(rho_a.*log(lA)); %technology shock process
e9=-log(Z)+(rho_z.*log(lZ)); %preference shock process
e10=-nu+(rho_nu.*lnu); %monetary policy shock process
e11=MC-(W_real./((1-alppha).*Y./N.*S)); %definition marginal cost
e12=1-(theta.*Pi.^(epsilon-1)+(1-theta).*(Pi_star).^(1-epsilon)); %LoM prices
e13=S-((1-theta).*Pi_star.^(-epsilon/(1-alppha))+theta*Pi.^(epsilon/(1-alppha)).*lS); %LoM price dispersion
e14=Pi_star.^(1+epsilon*(alppha/(1-alppha)))-(x_aux_1./x_aux_2*(1-tau)*epsilon/(epsilon-1)); %FOC price setting
e15=x_aux_1-(Z.*C.^(-siggma).*Y.*MC+betta*theta.*fPi.^(epsilon+alppha*epsilon/(1-alppha)).*fx_aux_1); %auxiliary price setting recursion
e16=x_aux_2-(Z.*C.^(-siggma).*Y+betta*theta*fPi.^(epsilon-1).*fx_aux_2); %auxiliary price setting recursion 2
e17=Pi-(P./lP); %definition inflation

%Put the equations together
system = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17];

% Useful vectors (v) not to confuse things. 
vXl = [lA lZ lS lPi_star lP lMC lR lrealinterest lPi lQ lN lC lW_real lY lx_aux_1 lx_aux_2 lnu]; %lagged variables
vX = [A Z S Pi_star P MC R realinterest Pi Q N C W_real Y x_aux_1 x_aux_2 nu]; %contemporaneous variables
vXf = [fA fZ fS fPi_star fP fMC fR frealinterest fPi fQ fN fC fW_real fY fx_aux_1 fx_aux_2 fnu]; %forward variables
vXss = repmat(xss',1,3)'; 
vXss = vXss(:)'; %steady state values for all variables
vVars = reshape([vXl;vX;vXf],size(vXl,1),[]); %all variables (all periods)

%% --------------------------Linearization--------------------------------
% Linearised system A x(t-1)+B x(t)+C x(t+1)=0

A = jacobian(system,vXl); %symbolic NumVarxNumVar
A = double(subs(A,vVars,vXss)); %double NumVarxNumVar, evaluated at SS
B = jacobian(system,vX); 
B = double(subs(B,vVars,vXss));
C = jacobian(system,vXf); 
C = double(subs(C,vVars,vXss));

% Convert to log-linear system
M = ones(NumVar,1)*xss;
A = A.*M;
B = B.*M;
C = C.*M;

%% --------------------------Time Iteration Solver--------------------------------
% Initialize
metric = 1;
F = 0;
% Run the algorithm until convergence
% Find  x_t=F x_{t-1}+Q u_t.
TimeStartLti=tic;
while metric>Tol
    F = -(B+C*F)\A;
    metric = max(max(abs([A+B*F+C*F*F])));
end
Q = -inv(B+C*F);
TimeEndLTI=toc(TimeStartLti);
disp(['Algorithm run time was ',num2str(TimeEndLTI),' seconds']); 

%% --------------------------IRFs--------------------------------
% Some specifications based on what shock we're looking at 
if strcmp(Shock,'TFP')
    ShockPos = 8; %note that this refers to the equation number that we're "shocking"
    ShockSize = sig_a; 
elseif strcmp(Shock,'Preference')
    ShockPos = 9;
    ShockSize = sig_z;
elseif strcmp(Shock,'MP')
    ShockSize = sig_nu;
    ShockPos = 10;
end
IrfName=strcat(Shock, ' Shock ');

% Generate the IRFs
vu = zeros(NumVar,1); 
vu(ShockPos) = 100*ShockSize; %ShockSize*abs(randn)
mirf(:,1) = Q*vu;
IrfHorAx = [0:IrfHor-1];
IrfPlotNum = (round((numel(VarPlot)-2)/2)*2+2)/2;

for t = 1:IrfHor-1
    mirf(:,t+1) = F*mirf(:,t); %iterate forward
end

% Plot the IRFs
figure('Name',IrfName)
for j=VarPlot
subplot(IrfPlotNum,2,find(VarPlot==j))
p1=plot(IrfHorAx,mirf(j,:),'-o','linewidth',1.5,'color',[0 0.45 0.75]);
set(gca,'FontSize',12,'fontname','times')
xlabel('Period','FontSize',12,'fontname','times')
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title(VarNames(VarPlot==j))
end

%% --------------------------Simulation--------------------------------
% This simulates a time path for the variables chosen (same as IRF)
% given a standard normal shock (to either MP, TFP or preferences)
% It's straightforward to amend the code to include multiple shocks
NumSimPeriods = 200;
mu = zeros(NumVar,NumSimPeriods);
mu(ShockPos,:) = randn(1,NumSimPeriods);
% u(8,:)=randn(1,NumSimPeriods);
% u(9,:)=randn(1,NumSimPeriods);
% u(10,:)=randn(1,NumSimPeriods);
vSim(:,1) = Q*mu(:,1);
for t = 1:NumSimPeriods-1
    vSim(:,t+1) = F*vSim(:,t)+Q*mu(:,t+1);
end

SimName = "Simulation based on"+" "+IrfName;

figure('Name',SimName)
for j=VarPlot
subplot(IrfPlotNum,2,find(VarPlot==j))
p1=plot(vSim(j,:),'linewidth',1.5,'color','k');
set(gca,'FontSize',12,'fontname','times')
xlabel('Period','FontSize',12,'fontname','times')
ylabel('% dev. from SS','FontSize',12,'fontname','times')
title(VarNames(VarPlot==j))
end

%% --------------------------Done------------------------------------
TimeEnd=toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']); 