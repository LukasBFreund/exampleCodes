%%=========================================================================
% Simple TANK model 
% Runs on Dynare 4.4.3 and higher
% C Cantore and LB Freund
%%=========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF ENDOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var 
    UCS           ${U^{C^S}}$                    (long_name='Marginal Utility of Consumption, Savers')
    UCH           ${U^{C^H}}$                    (long_name='Marginal Utility of Consumption, Hand-to-mouth')
    UHH           ${U^{H^H}}$                    (long_name='Marginal Utility of Leisure, Hand-to-mouth')
    CS            ${C^S}$                        (long_name='Consumption, Savers')
    CH            ${C^H}$                        (long_name='Consumption, Hand-to-mouth')
    HH            ${H^H}$                        (long_name='Hours, Hand-to-mouth')
    HS            ${H^S}$                        (long_name='Hours, Savers')
    UHS           ${U^{H^S}}$                    (long_name='Marginal Utility of Leisure, Savers')
    H             ${H}$                          (long_name='Aggregate Hours') 
    R             ${R}$                          (long_name='Real Interest Rate')
    Rn            ${R^n}$                        (long_name='Nominal Interest Rate')
    PIE           ${\Pi}$                        (long_name='Inflation')
    W             ${W}$                          (long_name='Real Wage')
    Y             ${Y}$                          (long_name='Real Output')
    MPL           ${MPL}$                        (long_name='Marginal Product of Labor')
    MC            ${MC}$                         (long_name='Real Marginal Costs')
    C             ${C}$                          (long_name='Consumption')
    Z             ${Z}$                          (long_name='Labor Augmenting shock process')
    profits       ${D}$                          (long_name='Profits - aggregate')
    profitsC      ${D^S}$                        (long_name='Profits - Savers')
    ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
    epsZ           ${\epsilon^{Z}}$       (long_name='Technology shock')
    epsM           ${\epsilon^{M}}$       (long_name='Monetary Policy shock')
    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters 
   betta          ${\beta}$               (long_name='Discount Factor') 
   sigma_c        ${\sigma}$              (long_name='Intertemporal elasticity of substitution')
   varrho         ${\varphi}$             (long_name='Inverse of Frish elasticity of Labor Supply')  
   alp            ${\alpha}$              (long_name='decreasing returns to scale')
   phi_p          ${\phi^p}$              (long_name='Rotemberg price adj. costs')
   zzeta          ${\epsilon}$            (long_name='Elasticity of substitutions between intermediate goods varieties')
   rhoZ           ${\rho^{Z}}$            (long_name='autoregressive parameter for Technology shock')
   rhoG           ${\rho^{G}}$            (long_name='autoregressive parameter for Government Spending shock')   
   theta_pie      ${\phi^{\pi}}$          (long_name='Taylor rule coeff of inflation')  
   gammap         ${\gamma^{p}}$          (long_name='Price Indexation') 
   lambda         ${\lambda}$             (long_name='Share of Hand-to-mouth Agents')
   tauD           ${\tau^D}$              (long_name='Tax on Profits')
   tauS           ${\tau^S}$              (long_name='Production Subsidy')
   psi            ${\psi}$                (long_name='Slope of Phillips Curve')
   //The following parameters are steady state realtionships
   Hss            ${\bar H}$              (long_name='Steady State Hours') 
   PIEss          ${\bar \Pi}$            (long_name='Steady State Inflation') 
   LSss           ${\bar{s}^{h}}$         (long_name='Steady State Labor Share') 
   ;

%%%%%%%%%%%%%%%%%%%%%%
%%PARAMETERS VALUES %%
%%%%%%%%%%%%%%%%%%%%%%
    betta      = 0.99; 
    rhoZ       = 0.75;  
    sigma_c    = 1;
    varrho     = 0.2; 
    gammap     = 0;
    lambda     = 0.5;
    theta_pie  = 1.5; 
    rhoG       = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%STEADY STATE RELATIONSHIPS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PIEss=1;
    Hss=1;
    zzeta=6;                %elasticity of substitution between differentiated goods 
    LSss = 0.67;            %Steady state Labor share if you want to allow for decreasing returns to scale
    alp = 0;                %1-LSss/MCss; %decreasing returns decreasing returns to scale
    tauD = 0;               %1-(1+varrho)/((1-lambda)^(-1)*varrho)-0.001;%Tax on Profits
    tauS = (zzeta-1)^(-1);  %Production Subsidy
    s_prices_duration = 3.5;
    calvo = 1-1/s_prices_duration;
    phi_p = calvo*(zzeta-1)/((1-calvo)*(1-betta*calvo)); % implied Rotemberg parameter (exploit first-order equivalence)
    psi = zzeta/phi_p; %slope phillips curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MODEL EQUATIONS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model(linear);

[name='Marginal Utility of Consumption, Savers']
UCS=-sigma_c*(CS);

[name='Euler Equation, Savers']
UCS=R+UCS(+1);

[name='Marginal Utility of Leisure, Savers']
UHS=varrho*HS;

[name='Labor Supply, Savers']
W=UHS-UCS;

[name=' Marginal Utility of Consumption, Hand-to-mouth']
UCH=-sigma_c*(CH);

[name='Marginal Utility of Leisure, Hand-to-mouth']
UHH=varrho*HH;

[name='Consumption, Hand-to-mouth']
CH=W+HH+tauD/lambda*profits;

[name='Labor Supply, Hand-to-mouth']
W=UHH-UCH;

[name='Cobb-Douglas Prodution function']
Y=(Z+H)*(1-alp);

[name='Marginal product of Labor']
MPL=H*(-alp)+(Z)*(1-alp);

[name='Real Wage']
W=MC+MPL;

[name='Profits']
profits=-MC;

[name='Profits - Savers']
profits=(1-lambda)*profitsC;

[name='Resource Constraint'] 
Y=C;

[name='Fisher Equation'] 
R=Rn-PIE(+1);

[name='Aggregation: Consumption']
C=lambda*CH+(1-lambda)*CS;

[name='Aggregation: Labor']
H=lambda*HH+(1-lambda)*HS;
%HH=HS;

[name='Linear Phillips Curve '] 
PIE=betta*PIE(+1)+psi*MC; 

[name='Taylor Rule'] 
Rn=theta_pie*PIE+epsM;

[name='Labor Augmenting Shock'] 
Z=rhoZ*Z(-1)+epsZ;

end;

%%%%%%%%%%%%%%%%%%%%%%
%%SHOCKS            %%
%%%%%%%%%%%%%%%%%%%%%% 

shocks;
var epsZ; stderr 1;
var epsM; stderr 1;
end;

steady;
check;
resid(1);

/*
write_latex_parameter_table;
write_latex_dynamic_model;
write_latex_definitions;
collect_latex_files;
*/

%%%%%%%%%%%%%%%%%%%%%%%%%
%%STOCHASTIC SIMULATION%%
%%%%%%%%%%%%%%%%%%%%%%%%%
stoch_simul(order=1,irf=20,nograph,noprint);

