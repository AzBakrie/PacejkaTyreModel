clear
close all
clc

addpath(genpath('TTCdata'))

tyre_model = struct;
tyre_model.p = struct; % all the parameters will be stored here 

% load all the different runs 
run1 = load('B1965raw12.mat'); % latitude data
run2 = load('B1965raw13.mat'); % latitude data

run3 = load('B1654run35.mat'); % longitude data
run4 = load('B1654run36.mat'); % longitude data


%% FY lateral slip
tyre_model.fy0 = FYclass(run1, run2);
tyre_model.fy0 = tyre_model.fy0.createDataBins(false);
tyre_model.fy0 = tyre_model.fy0.create_sweeps(1000, false);
tyre_model.fy0 = tyre_model.fy0.fitCubicSpline();
tyre_model.fy0 = tyre_model.fy0.generateSplinePoints();
[tyre_model.fy0, tyre_model.p] = tyre_model.fy0.magicFormulaParamterOptimisation(tyre_model.p);
tyre_model.fy0 = tyre_model.fy0.buildPlotLibrary("combined", "no");
tyre_model.fy0 = tyre_model.fy0.resultsAnalysis(false);

%% FX longitude slip
tyre_model.fx0 = FXclass(run3, run4);
tyre_model.fx0 = tyre_model.fx0.createDataBins();
tyre_model.fx0 = tyre_model.fx0.create_sweeps(500, false);
tyre_model.fx0 = tyre_model.fx0.fitCubicSpline();
tyre_model.fx0 = tyre_model.fx0.generateSplinePoints();
[tyre_model.fx0, tyre_model.p] = tyre_model.fx0.magicFormulaParamterOptimisation(tyre_model.p);
tyre_model.fx0 = tyre_model.fx0.buildPlotLibrary("combined", "no");
tyre_model.fx0 = tyre_model.fx0.resultsAnalysis(false);

%% MZ lateral slip
tyre_model.mz0 = MZclass(run1, run2);
tyre_model.mz0 = tyre_model.mz0.createDataBins();
tyre_model.mz0 = tyre_model.mz0.create_sweeps(1000, false);
tyre_model.mz0 = tyre_model.mz0.fitCubicSpline();
tyre_model.mz0 = tyre_model.mz0.generateSplinePoints();
[tyre_model.mz0, tyre_model.p] = tyre_model.mz0.magicFormulaParamterOptimisation(tyre_model.p);
tyre_model.mz0 = tyre_model.mz0.buildPlotLibrary("combined", "no");
tyre_model.mz0 = tyre_model.mz0.resultsAnalysis(false);

%% FY combined slip
tyre_model.fy = FYclass(run3, run4);
tyre_model.fy = tyre_model.fy.createDataBins(true);
tyre_model.fy = tyre_model.fy.create_sweeps(300, true);
tyre_model.fy = tyre_model.fy.fitCombinedCubicSpline();
tyre_model.fy = tyre_model.fy.generateCombinedSplinePoints(); 
[tyre_model.fy, tyre_model.p] = tyre_model.fy.magicFormulaCombinedParameterOptimisation(tyre_model.p);
tyre_model.fy = tyre_model.fy.buildPlotLibrary("combined", "yes");
tyre_model.fy = tyre_model.fy.resultsAnalysis(true);

%% FX combined slip
tyre_model.fx = FXclass(run3,run4);
tyre_model.fx = tyre_model.fx.createDataBins(true);
tyre_model.fx = tyre_model.fx.create_sweeps(300, true);
tyre_model.fx = tyre_model.fx.fitCombinedCubicSpline();
tyre_model.fx = tyre_model.fx.generateCombinedSplinePoints();
[tyre_model.fx, tyre_model.p] = tyre_model.fx.magicFormulaCombinedParameterOptimisation(tyre_model.p);
tyre_model.fx = tyre_model.fx.buildPlotLibrary("combined", "yes");
tyre_model.fx = tyre_model.fx.resultsAnalysis(true);

%% MZ COMBINED SLIP - NEED TO FILTER MZ VALS
tyre_model.mz = MZclass(run3, run4);
tyre_model.mz = tyre_model.mz.createDataBins(true);
tyre_model.mz = tyre_model.mz.create_sweeps(300,true);
tyre_model.mz = tyre_model.mz.fitCombinedCubicSpline();
tyre_model.mz = tyre_model.mz.generateCombinedSplinePoints();
tyre_model.mz.evaluation_points.all_Fy = tyre_model.fy.magic_formula.F_fit ;
tyre_model.mz.evaluation_points.all_Fx = tyre_model.fx.magic_formula.F_fit ; 
[tyre_model.mz, tyre_model.p] = tyre_model.mz.magicFormulaCombinedParameterOptimisation(tyre_model.p);
tyre_model.mz = tyre_model.mz.buildPlotLibrary("combined", "yes");
tyre_model.mz = tyre_model.mz.resultsAnalysis(true);

%% SAVING 
p = tyre_model.p ;
save('p.mat', 'p')

% writetable(struct2table(tyre_model.p), 'p.csv')
