XAS_analysis: SES1 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess, mmap, numpy(version 1.13.1). It was developed on a windows system using Anaconda.


FILES INCLUDED 

1,) README.md				This currently holds all the documentation for
							this repository

2,) SES1.py  				Python source code for the excited state (theoretical)
							calculations of a liquid
						
3,) Imidazole_geom.xyz		A molecular geometry file used for testing 
							and development. The name of this file is used to identify 
							results directories so it is useful to give this a name that 
							is relevant to the chemicals being analysed.
						
4,) orca_parameters.txt 	This text file contains an orca keyword for single point 
							energy calculations on each line. It allows the user to 
							overwrite the default ocra values that are in this program.

							
INFORMATION ON THE LES1 (Liquid THEORETICAL EXCITED STATE CALCULATION) CODE

usage: 

LES1.py [-h] [-geom_file_name FILE] [-opi FILE] [-orca] [-path_out] [-element] [-h_opt] DIR


GTE1: Theoretical electron density function calulation

positional arguments:

  DIR					The Directory containing molecular geometry files to be read in; 
						Directory name must include its path and it is best not to be in 
						the directory structure for this code.

optional arguments:

  -h, --help          	show this help message and exit
  
  -geom_file_name FILE  If only one geometry file is used then molecular geometry file 
						name should be provided rahter than direcoty. File name must 
						include its path and it is best not to be in the directory 
						structure for this code
					  
  -opi FILE           	input file with orca optimation parameters. This over
						write any default orca optimisation parameters.
					  
  -orca 	          	full path to orca program; C:\Orca\orca is the
						default path
					
  -path_out				The directory where the results directory will be 
						written

  -element				The elemtent for which excited state calculations 
						will be performed
	
  -h_opt				Determines whether to run geometry optimization calculation 
						for hydrogen positions 

EXAMPLE USAGE

1, minimal needed to run:

GTE1.py C:\Users\userid\Desktop\data\Imidazole\Imidazole_geom.xyz

2, minimal needed to run with user defined orca parameters:

GTE1.py C:\Users\userid\Desktop\data\Imidazole\Imidazole_geom.xyz -opi C:\Users\userid\Desktop\analysis\Imidazole\orca_parameters.txt


BACKGROUND

This program was developed To work with the software ORCA version 4.0.1.2 which can be downloaded after you have registered at the ORCA forum:

https://orcaforum.cec.mpg.de/


This program is the SES1 part of the computational pipeline or workflow. This is the Theory Excited State calculation pipeline for solids. The computation is in 2 steps:
1, Single Point Energy Calculation or Hydrogen Position Geometry Optimization 2, Excited stated calculation or TDDFT calculation (Time Dependent Density Function Theory)

If the software finishes early and not all of the steps complete sucessfully there will be an error message printed. These error messages are:
ERROR 1: ORCA did not terminate normally (to be implemented)
