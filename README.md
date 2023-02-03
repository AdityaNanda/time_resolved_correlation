# Time resolved correlation of distributed brain activity tracks E-I balance and accounts for diverse scale-free phenomena

This is a MATLAB based tool for sampling timeseries data with preserved mean, variance and time-resolved correlation. Time-resolved correlation, defined as Pearson's correlation between adjacent timepoints, has been shown to track and recapitulate 1/f-based estimates of global cortical Excitation-Inhibition (E-I) balance. The tool is based on nullspace sampling methods. Although, primarily designeed with electrophysiologgy datasets in mind, the tool can be used for any timeseries data. 

# Getting started

Download or clone this repository into your preferred location. To download, simply press the _Code_ button and select _Download ZIP_.  To clone, enter `git clone https://github.com/AdityaNanda/time_resolved_correlation` on the command line.

To add the toolbox to the MATLAB path, use the command `addpath(genpath(path_to_toolbox))` or use the _Set path_ button in the _Environment_ section of the _Home_ ribbon, and click _add with subfolders_. Now you can directly access the relevant functions.

# Demo

The script 'main.m' includes demo of the two principal functions and uses ~ 20 seconds of freely available Human iEEG data [1] to demostrate their use . The two functions are 

*  'mean_var_corr.m'  - > generate synthetic timeseries with preserved time-resolved correlation of distributed brain activity (and mean and varaince). 
* 'mean_var.m'        - > generate synthetic timeseries with preserved mean and variance

# Reference

[1] Elizabeth Johnson (2018); Intracranial EEG recordings of medial temporal, lateral frontal, and orbitofrontal regions in 10 human adults performing a visuospatial working memory task. CRCNS.org <a href="http://dx.doi.org/10.6080/K0VX0DQD">doi</a>

