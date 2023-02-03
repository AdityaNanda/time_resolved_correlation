# Time resolved correlation of distributed brain activity tracks E-I balance and accounts for diverse scale-free phenomena

This is a set of MATLAB-based tools for sampling timeseries data with preserved mean, variance and time-resolved correlation. Time-resolved correlation is the Pearson correlation between brain activity vectors at adjacent timepoints and tracks 1/f-based estimates of excitation-inhibition balance. The present tool leverages nullspace sampling methods. Although it is primarily designeed with electrophysiology datasets in mind, it can be used for any multivariate timeseries data. 

These tools also contain functions for computing detrended fluctuation analysis and avalanche statistics.

# Getting started

Download or clone this repository into your preferred location. To download, simply press the _Code_ button and select _Download ZIP_.  To clone, enter `git clone https://github.com/AdityaNanda/time_resolved_correlation` on the command line.

To add the toolbox to the MATLAB path, use the command `addpath(genpath(path_to_toolbox))` or use the _Set path_ button in the _Environment_ section of the _Home_ ribbon, and click _add with subfolders_. Now you can directly access the relevant functions.

# Demo

The script `demo_main.m` is a demo of the two main functions on a ~20 second recording of a freely available human intracranial EEG dataset [(Johnson, 2018)](http://dx.doi.org/10.6080/K0VX0DQD). The two functions are:

* `mean_var_corr.m`: sample timeseries constrained by time-resolved correlation, variance, and mean. 
* `mean_var.m`: sample timeseries constrained by time-resolved variance and mean.
