## Analysis of Rodent Respiration
### Objective
The objective of this project is to analyze the rodent respriation signals from medullary neurons and actual breathing rhythms by detecting their maxima and minima.
### Description
#### Loading .abf files
The '.abf' files that include raw and pre-processed EEG and breathing signals can be loaded into MATLAB using 'abfload.m'. This code was made by Jo Suresh and Tahra Eissa in 2016.
#### Peak Detection and Averaging
The 'Rodent_breathing_analysis.m' file filters the signal in both low-pass and high-pass and detects maxima and minima of the selected signals. This file is identical to the 'Rodent_breathing_analysis_w_avg.m' file, except that the latter adds a short extension, which enables amplitude averaging of detected peaks over the time window defined by the user. The 'avg_peaks.m' function is employed for peak averaging, and the 'manualedit_maxmin.m' allows user to manually edit the initial output data by clicking on the generated plots.
#### Requirements
MATLAB 2018 b
.abf data
### Acknowledgments
The set of codes above was written based on previously developed code of Dr. Eissa and Suresh and by consulting Prof. Tryba at the University of Chicago Department of Pediatrics Section of Neurology.
