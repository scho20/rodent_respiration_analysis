## Analysis of Rodent Respiration Data
### Objective
The objective of this project is to analyze the rodent respriation signals from medullary neurons and actual breathing rhythms by detecting their maxima and minima.
### Description
#### 1. Loading .abf files
The '.abf' files that include raw and pre-processed EEG and breathing signals can be loaded into MATLAB using '[abfload.m](https://https://github.com/scho20/rodent_respiration_analysis/blob/master/abfload.m)'. This code was made by Jo Suresh and Tahra Eissa in 2016.
#### 2. Peak Detection and Averaging
The '[Rodent_breathing_analysis.m](https://https://github.com/scho20/rodent_respiration_analysis/blob/master/Rodent_breathing_analysis.m)' file filters the signal in both low-pass and high-pass and detects maxima and minima of the selected signals. This file is identical to the '[Rodent_breathing_analysis_w_avg.m](https://github.com/scho20/rodent_respiration_analysis/blob/master/Rodent_breathing_analysis_w_avg.m)' file, except that the latter adds a short extension, which enables amplitude averaging of detected peaks over the time window defined by the user. The '[avg_peaks.m](https://github.com/scho20/rodent_respiration_analysis/blob/master/avg_peaks.m)' function is employed for peak averaging, and the '[manualedit_maxmin.m](https://github.com/scho20/rodent_respiration_analysis/blob/master/manualedit_maxmin.m)' allows user to manually edit the initial output data by clicking on the generated plots.
#### 3. Requirements
MATLAB 2018b <br/>
'.abf' data files
### Acknowledgments
The codes above were written based on previously developed ones of Dr. Eissa and Suresh and by consulting Prof. Tryba at the University of Chicago Department of Pediatrics.
