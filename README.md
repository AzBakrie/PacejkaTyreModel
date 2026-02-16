# PacejkaTyreModel
Pacejka Tyre Model on MATLAB for Imperial Formula Racing 

MATLAB framework to process the TTC tyre test data and identify the parameters for Pacejka Magic Formula.
The toolkit filters raw tyre test runs, extracts the steady state sweeps, fits smoothing splines and optimises the Magic Formula coefficients for:
- Lateral Force (Fy)
- Longitudinal Force (Fx)
- Aligning Moment (Mz)
- Combined Slip Behaviour

# Pipeline Overview 
Raw TTC Data\
  → Condition Binning \
  → Sweep Detection\
  → Spline Fitting\
  → Evaluation Point Generation\
  → Parameter Optimisation\
  → Magic Formula Model\
  → Plotting\

# Requirements
MATLAB toolbox requires:
- Optimization Toolbox
- Curve Fitting Toolbox (csaps)

# Installation 
Clone the repository\
git clone https://github.com/yourusername/tyre-model.git cd tyre-model \
Open MATLAB and add to path:\
addpath(genpath(pwd)) 

# Basic Usage
Run the mainscript.m file

# Plotting
After running the mainscript.m file to show the plots, write in the Command Window: \
tyre_model.obj.showPlot(load_value)\
Where load value is the load you would like to display 






