XAS_analysis (X-ray Absorption Spectroscopy)

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1) and timeit. It was developed on a windows system using Anaconda.

SOFTWARE OVERVIEW

This software is the start of automating the computational pipeline or workflow for XAS Analysis. There are python 4 scripts for each state of matter ie, gas, liquid and solid, however the experimental spectra peak fitting (E2) and the comparison (C1) are the same for all these states. The initial of each of these states is used at the start of a script name to indicate what the script is used for. The scripts that make up this pipeline is in its own directory. These are:

E1: Experimental spectra background noise subtraction and normalization. This is still to be implemented as a deep learning problem.

E2: Experimental spectra peak fitting

GTE1: Theoretical electron density function calulation for a gas

LES1: (Theoretical) excited state calulation for a liquid

SES1: (Theoretical) excited state calulation for a solid

C1:  Compares experimental spectra with theorectically calculated data

GW1:  Wrapper Script to run E2, GTE1 and C1 for a gas

LW1:  Wrapper Script to run E2, LES1 and C1 for a liquid

SW1:  Wrapper Script to run E2, SES1 and C1 for a solid


The data file edge_data.txt which is a text file containing edge data information (derived from "Journal of Physics: Conference Series 712 (2016) 012070" which is taken from the X-ray Data Booklet) is also stored in this directory as it is used by more than one script.

NOTES TO DEVELOPERS
There are 2 technical choices made during development that need explaining because of the limitations of varies systems at the time the software was developed.

Firstly the issue of the html reports that are created automatically when some scripts run. Images in these scripts are created and are referenced in an html report. The images exist separately to this report and can be used in other documents. The resolution of these images is determined by the graphics settings of the system that creates them. These are scaled in the html to 100% so that high resolution images will be visible on the screen (it will also increase the size of low resolution images). At the time of writing the html 5 format was new and was not fully implemented in all browsers. Also the html 5 standard had new ways to handle images that did not include scaling so the html template is not compatible with html 5.

Secondly, all the scripts create log files that provide meta data on when each sctipt is executed. These are plain text files and do not use the log file handler. This is because at the time of development the arparse library was not compatible with the log file handler.








