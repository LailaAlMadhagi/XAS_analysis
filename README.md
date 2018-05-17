#XAS_analysis 
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17 The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1) and timeit

This program was developed with work with the software ORCA version 4.0.1.2 which can be downloaded after you have registered at the ORCA forum:
https://orcaforum.cec.mpg.de/

This program is the TE1 part of the computational pipeline or workflow. This is the Theory Electronic (Structutre) 1 pipe. The computation is in 3 steps:
1, Geometry Optimisation 2, Frequency Calculation 3, Excited stated calculation or TDDFT calculation (Time Dependent Density Function Theory)

The sucess of the execution of this software relies on each step sucessfully completing. The first step, the optimisation step) may require a variety of optimisation strategies to be applied before it is successful. The software applies a list of strategies until success is reached or until there are no more suitable strategies available. Each of these strategies relates to a setting in ORCA. These can be set in a pre-defined order within the code or the user can create their own list in an input file (this is in the process of being implemented).

If the software finishes early and not all of the steps complete sucessfully there will be an error message printed. These error messages are:
ERROR 1: ORCA did not terminate normally
