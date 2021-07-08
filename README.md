# noisy-cell-growth

**Paper Title**:
Cell growth model with stochastic gene expression helps understand the growth advantage of metabolic exchange and auxotrophy

**Authors**:
Dibyendu Dutta #, Supreet Saini

Department of Chemical Engineering, Indian Institute of Technology Bombay, Mumbai, India

**Running Title**: How do cooperators grow faster than prototrophs 

**# Corresponding Authors**: 
Dibyendu Dutta	Email: dibyendu_dutta@iitb.ac.in

==============================================================================================

Simulations carried out in MATLAB 2019b.

The folder "modified_added_sizer", must be first added to the path in MATLAB, using code below, since it contains many of the functions that are used for the simulation runs.
addpath('path_to_\modified_added_sizer')

MATLAB Simulation code used in the manuscript are in the folder Manuscript Simulations. The subfolders are named based on the figures.
In each folder there are files names such as "run_xxx.m"  which can be opened and edited to change the run parameters and the location of simulation output.

The file "analysis_xxx.m", or or "A_xxx.m" contain the necessary codes for analysing the data generated from the previous steps. These files need to be edited to give the location of the generated data files, and accordingly the part of the code required must be executed. The analysis code is not a one step run script, but a collection of short scripts.

For any doubts please reach out to the email given above.

Best,
Dibyendu Dutta
