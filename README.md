# PacejkaTyreModel
Pacejka Tyre Model on MATLAB for [Imperial Formula Racing](https://www.imperialformularacing.com/)

MATLAB framework to process the TTC tyre test data and identify the parameters for Pacejka Magic Formula.
The toolkit filters raw tyre test runs, extracts the steady state sweeps, fits smoothing splines and optimises the Magic Formula coefficients for:
- Lateral Force (Fy)
- Longitudinal Force (Fx)
- Aligning Moment (Mz)
- Combined Slip Behaviour

This project was made possible with the data provided by [Formula SAE Tire Test Consortium (FSAE TTC)](https://www.millikenresearch.com/fsaettc.html) and [Calspan Tire Testing Facility (TIRF)](https://calspan.com/automotive/fsae-ttc)

# Prior Reading
It is recommended that users read the documents in the [TTCRunGuide](https://github.com/AzBakrie/PacejkaTyreModel/tree/main/TTCRunGuide) Folder and to register in the [TTC Private Forum](https://www.fsaettc.org/)

# Pipeline Overview 
Raw TTC Data\
  → Condition Binning \
  → Sweep Detection\
  → Spline Fitting\
  → Evaluation Point Generation\
  → Parameter Optimisation\
  → Magic Formula Model\
  → Plotting

# Requirements
MATLAB toolbox requires:
- Optimization Toolbox
- Curve Fitting Toolbox (csaps)

# Installation 
Clone the repository
```bash
git clone https://github.com/yourusername/tyre-model.git
cd tyre-model
```

Open MATLAB and add to path:
```matlab
addpath(genpath(pwd)) 
```
# Basic Usage
Run the ```mainscript.m``` file, ensuring the other classes are in the same path 

  # Plotting
  After running the mainscript.m file to show the plots, write in the Command Window:
  ```matlab
  tyre_model.obj.showPlot(load_value)
  ```
  Where ```load_value``` is the load you would like to display 






