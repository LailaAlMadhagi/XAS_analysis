XAS_analysis (X-ray Absorption Spectroscopy)

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1) and timeit. It was developed on a windows system using Anaconda.

SOFTWARE OVERVIEW

This software is the start of automating the computational pipeline or workflow for XAS Analysis. There are 6 python 3 scripts that make up this pipeline and each of these is in its own directory. These are:

E1: Experimental spectra background noise subtraction and normalization

E2: Experimental spectra peak fitting

TE1: Theoretical electron density function calulation

C1:  Compares experimental spectra with theorectically calculated data.

TMD1: Creates a structural model of the system

TES1: A Theortical Electronic Structure Excited State Calculation for a Solution


The data file edge_data.txt which is a text file containing edge data information (derived from "Journal of Physics: Conference Series 712 (2016) 012070" which is taken from the X-ray Data Booklet) is also stored in this directory as it is used by more than one script.









