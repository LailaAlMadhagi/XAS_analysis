XAS_analysis: GTE1 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1). It was developed on a windows system using Anaconda.


FILES INCLUDED 

1,)README.md			This currently holds all the documentation for
						this repository

2,) GTE1.py  			Python source code for the theoretical
						electronic structure calculations of a gas
						
3,) Imidazole_geom.xyz	A molecular geometry file used for testing 
						and development. The name of this file is used to identify results directories so it is useful to give this a name that is relevant to the chemicals being analysed.
						
4,) orca_parameters.txt This text file contains an orca keyword for
						optimisation on each line. It allows the user to overwrite the default ocra values that are in this program.


						
INFORMATION ON THE GTE1 (GAS THEORETICAL ELECTRONIC STRUCTURE CALCULATION) CODE

usage: 

GTE1.py [-h] [-opi FILE] [-orca] [-path_out] [-element] FILE


GTE1: Theoretical electron density function calulation

positional arguments:

  FILE                Molecular geometry file to be read in; filename must include
                      its path and it is best not to be in the directory structure for this code.

optional arguments:

  -h, --help          	show this help message and exit
  
					  
  -opi FILE           	input file with orca optimation parameters. This over
						write any default orca optimisation parameters.
					  
  -orca 	          	full path to orca program; C:\Orca\orca is the
						default path
					
  -path_out				The directory where the results directory will be 
						written

  -element				The elemtent for which excited state calculations 
						will be performed
	

EXAMPLE USAGE

1, minimal needed to run (defaults assume it is a gas):

GTE1.py C:\Users\userid\Desktop\data\Imidazole\Imidazole_geom.xyz

2, minimal needed to run with user defined orca parameters:

GTE1.py C:\Users\userid\Desktop\data\Imidazole\Imidazole_geom.xyz -opi C:\Users\userid\Desktop\analysis\Imidazole\orca_parameters.txt


BACKGROUND

This program was developed To work with the software ORCA version 4.0.1.2 which can be downloaded after you have registered at the ORCA forum:

https://orcaforum.cec.mpg.de/


This program is the GTE1 part of the computational pipeline or workflow. This is the Theory Electronic (Structutre) 1 pipe for gas. The computation is in 3 steps:
1, Geometry Optimisation 2, Frequency Calculation 3, Excited stated calculation or TDDFT calculation (Time Dependent Density Function Theory)

The sucess of the execution of this software relies on each step sucessfully completing. The first step, the optimisation step) may require a variety of optimisation strategies to be applied before it is successful. The software applies a list of strategies until success is reached or until there are no more suitable strategies available. Each of these strategies relates to a setting in ORCA. These can be set in a pre-defined order within the code or the user can create their own list in an input file (this is in the process of being implemented).

If the software finishes early and not all of the steps complete sucessfully there will be an error message printed. These error messages are:
ERROR 1: ORCA did not terminate normally
