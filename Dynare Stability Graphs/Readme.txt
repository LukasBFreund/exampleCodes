%% Readme
%--------------------------------------------------------------------------
This code creates stability maps based on the Dynare toolbox, i.e., regions of saddle path stability, indeterminacy and instability.
Runs on Dynare 4.4.3 and higher
LB Freund (luk.freund@gmail.com)
Last updated: November 2019

%% Content
%--------------------------------------------------------------------------
Files included:
i.) PlotStabilityGraphs.m: main file.
ii.) TANKs_linear.mod: example Dynare model code called by i.).

%% Example run
%--------------------------------------------------------------------------
1. As "ModelName" choose the name of the Dynare file (w/o ".mod").
2. Pick two parameters as "Parameter1" and "Parameter2". These need to have the same name as the corresponding parameter values in the Dynare code.
3. Choose min/max/steps to define the space over which the parameters are evalued.
4. Optionally: enable "OptionPrint" and choose a file name and output path.
5. Run.


%% Additional notes
%--------------------------------------------------------------------------
- The key steps in PlotStabilityGraphs.m are lines 54-76.
- TANKs_linear is a barebones two-agent New Keynesian model (only labor).
- The rounding of the axes labels is just a preference and can be adjusted in lines 126-130 & 135-140.

