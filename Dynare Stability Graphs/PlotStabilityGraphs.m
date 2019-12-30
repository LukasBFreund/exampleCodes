%%=========================================================================
% This code creates stability maps based on the Dynare toolbox, i.e.,
% regions of saddle path stability, indeterminacy and instability.
% Runs on Dynare 4.4.3 and higher
% LB Freund
% Last updated: November 2019
%%=========================================================================
 
%% Information
%--------------------------------------------------------------------------
% This approach uses the following output from Dynare's built-in resol function:
%   info(1) = 0     =>    No error.
%   info(1) = 3     =>    Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%   info(1) = 4     =>    Blanchard & Kahn conditions are not satisfied: indeterminacy.

%% Housekeeping
%--------------------------------------------------------------------------
clear;
close all;
clc;
TimeStart = tic;

%% User Choices
%--------------------------------------------------------------------------
% Choose model (name of .mod file)
ModelName = 'TANKs_linear';

% Set parameter space to be evaluated
Parameter1 = 'varrho'; 
Parameter1Min = 0.01;
Parameter1Max = 20;
Parameter1Steps = 50;

Parameter2 = 'lambda';
Parameter2Min = 0.01;
Parameter2Max = 1;
Parameter2Steps = 50;

% Axes labels
Parameter1String = horzcat('\',Parameter1); % assumes Parameter1 is written consistent w/ latex notation, else choose manually
Parameter2String = horzcat('\',Parameter2); 

% Select printing options
OptionPrint = 1;        
FigName = 'Fig_Stability_Spec';
TargetPath = '.\Output\';   

%% Run Dynare
%--------------------------------------------------------------------------
%dynare TANKs_linear noclearall 
Model=['dynare ' ModelName ' noclearall'];
eval(Model)

%% Compute Stability Regions
%--------------------------------------------------------------------------
% Create parameter vectors
vParameter1 = linspace(Parameter1Min,Parameter1Max,Parameter1Steps);
vParameter2 = linspace(Parameter2Min,Parameter2Max,Parameter2Steps);
[mParameter1,mParameter2] = meshgrid(vParameter1,vParameter2); % matrix version

% Initialize
Iter = 0;
options_.qz_criterium = 1+1e-6;
mStability = NaN(size(mParameter1));

% Compute stability regions
for iP1 = 1:length(vParameter1)
    for iP2 = 1:length(vParameter2)
       Iter = Iter+1;
       set_param_value(Parameter1,mParameter1(iP1,iP2)); 
       set_param_value(Parameter2,mParameter2(iP1,iP2)); 
       [dr,info] = resol(0,M_,options_,oo_); % compute perturbation based decision rules 
       mStability(iP1,iP2) = info(1);        % store stability info
    end
end

%% Plot
%--------------------------------------------------------------------------
%% Design defaults 
OptionGreycolor = 0;    % if want grey colorscheme instead of standard colors
vLinestyle = {'-','-.','--',':'};
FontsizeDefault = 10;
FontsizeAxis = 10;
FontSizeLegend = 8;
FontsizeAxisticks = 8;
Fonttype = 'times';
LinewidthDefault = 1.6;
LinewidthAlt = 1;
ColorZeros = 'k';
StyleZeros = '-';
if OptionGreycolor == 0
vColors = {[4,30,150]/255,[1 0.5 0],[152,58,68]/255,[30,144,250]/255};
elseif OptionGreycolor == 1
vColors = {[0.2,0.2,0.2],[0.46,0.46,0.46],[0.6,0.6,0.6],[0.7 0.7 0.7]};
end

%% Figure 1: Simple contour plot focusing just on stability vs. not (not printed, just for illustration)
Fig1 = figure(1);
mDetplot = zeros(size(mStability));
mDetplot(mStability==0)=1; % cases when there's no error
figure(1)
contour(mParameter1,mParameter2,mDetplot,1,'k','Linewidth',LinewidthDefault)
xlabel(horzcat('$',Parameter1String,'$'),'fontsize',FontsizeAxis,'interpreter','latex'); 
ylabel(horzcat('$',Parameter2String,'$'),'fontsize',FontsizeAxis,'interpreter','latex');    
title('Determinacy Region','fontname',Fonttype,'Fontsize',FontsizeDefault);

%% Figure 2: Mesh plot (useful esp. if have all three cases)
TickScaleFactor = 10;       % factor reducing frequency of x-ticks relative to # steps 
Parameter1Stepsize = vParameter1(2)-vParameter1(1);
Parameter2Stepsize = vParameter2(2)-vParameter2(1);

Fig2 =figure(2);
markersize=1;
spy(mStability(:,:)==0,'r.');
set(get(gca,'children'),'color',vColors{2})
hold on 
spy(mStability(:,:)==4,'k.');
hold on 
spy(mStability(:,:)==3,'y.');
set(get(gca,'children'),'markersize',10)
axis xy;

% But need to scale axes labels to reflect parameter values rather than vec length
temp = {Parameter2Min:TickScaleFactor*Parameter2Stepsize:Parameter2Max};
% The following lines are just to make the graphs prettier, may wish to adjust
% depending on preferences
if Parameter2Max>10
N = 0;
else
N = 1;
end
temp = cellfun(@(x)round(x,N),temp,'UniformOutput',false);  % round
set(gca,'YTick',1:TickScaleFactor:length(vParameter2),'fontsize',FontsizeAxis); % this is important to show the right para values
set(gca,'YTicklabel',temp);

temp = {Parameter1Min:TickScaleFactor*Parameter1Stepsize:Parameter1Max};
if Parameter1Max>10
N = 0;
elseif Parameter1Max
N = 1;
end
temp = cellfun(@(x)round(x,N),temp,'UniformOutput',false);
set(gca,'XTick',1:TickScaleFactor:length(vParameter1),'fontsize',FontsizeAxis);
set(gca,'XTicklabel',temp);
ylabel(horzcat('$',Parameter2String,'$'),'fontsize',FontsizeAxis,'interpreter','latex');    
xlabel(horzcat('$',Parameter1String,'$'),'fontsize',FontsizeAxis,'interpreter','latex'); 

% Legend
legend1 = legend('Determinacy','Indeterminacy','Instability');
set(legend1,'fontname','times','Location','best','FontSize',FontSizeLegend,'Orientation','horizontal');
newPosition = [0.2316    0.9474    0.5518    0.0417];
newUnits = 'normalized';
set(legend1,'Position', newPosition,'Units', newUnits,'Orientation','horizontal');


%% Print
%--------------------------------------------------------------------------
if OptionPrint == 1
xSize = 17.5/2; 
ySize = 10; 
xCut = 0;
yCut = -0.5;

 set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')
  FigNamepdf =horzcat(horzcat(TargetPath,FigName),'.pdf');
print(FigNamepdf,'-dpdf','-painters')
end

%% Done!
%----------------------------------------------------------------------------
TimeEnd = toc(TimeStart);
disp(['Total run time was ',num2str(TimeEnd),' seconds']);