XAS_analysis: TE1 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1) and timeit. It was developed on a windows system using Anaconda.


FILES INCLUDED 

1,)README.md			This currently holds all the documentation for
						this repository

2,) TE1.py  			Python source code for the theoretical
						electronic structure calculations
						
3,) edge_data.txt		A text file containing edge data information
						(derived from "Journal of Physics: Conference Series 712 (2016) 012070" which is taken from the X-ray Data Booklet)
						
4,) geom.xyz            A molecular geometry file used for testing 
						and development. The name of this file is used to identify results directories so it is useful to give this a name that is relevant to the chemicals being analysed.
						
5,) orca_parameters.txt This text file contains an orca keyword for
						optimisation on each line. It allows the user to overwrite the default ocra values that are in this program.
						
6,) workflow.jpg 		An image file containing a diagram of all the
						planned pipes in this workflow.


						
INFORMATION ON THE TE1 (THEORETICAL ELECTRONIC STRUCTURE CALCULATION) CODE

usage: 

TE1.py [-h] [-op {gas,solution}] [-opi FILE] [-orca FILE] FILE


TE1: Theoretical electron density function calulation

positional arguments:

  FILE                molecular geometry file to be read in

optional arguments:

  -h, --help          show this help message and exit
  
  -op {gas,solution}  Select the default set of orca parameters for particular
                      chemical states.
					  
  -opi FILE           input file with orca optimation parameters. This over
                      write any default orca optimisation parameters.
					  
  -orca FILE          path to the orca executable; C:\Orca\orca.exe is the
                      default path


This program was developed TO work with the software ORCA version 4.0.1.2 which can be downloaded after you have registered at the ORCA forum:

https://orcaforum.cec.mpg.de/


This program is the TE1 part of the computational pipeline or workflow, a diagram of all the planned pipes is provided in the respository. This is the Theory Electronic (Structutre) 1 pipe. The computation is in 3 steps:
1, Geometry Optimisation 2, Frequency Calculation 3, Excited stated calculation or TDDFT calculation (Time Dependent Density Function Theory)

The sucess of the execution of this software relies on each step sucessfully completing. The first step, the optimisation step) may require a variety of optimisation strategies to be applied before it is successful. The software applies a list of strategies until success is reached or until there are no more suitable strategies available. Each of these strategies relates to a setting in ORCA. These can be set in a pre-defined order within the code or the user can create their own list in an input file (this is in the process of being implemented).

If the software finishes early and not all of the steps complete sucessfully there will be an error message printed. These error messages are:
ERROR 1: ORCA did not terminate normally
