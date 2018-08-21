XAS_analysis: E2 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/23. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later. 

This program was developed with Python 3.5.2 with modules: numpy(version 1.13.1), scipy(version 0.19.1), matplotlib(version 2.0.2) and lmfit(version 0.9.10).It was developed on a windows system using Anaconda

to install lmfit, follow the instructions on this web page:

https://anaconda.org/conda-forge/lmfit 


FILES INCLUDED

1,)E2.py						Python source code for peak fitting 
								
2,)edge_data.txt 				A text file containing edge data information
								(derived from "Journal of Physics: Conference Series 712 (2016) 012070" which is taken from the X-ray Data Booklet) 
								
3,) N1s_Imidazole_ISEELS.txt 	A text file containing experimental nitrogen K-edge data of imidazole in the gas phase. The file is used for testing and development 
								(data published in Apen et al, J. Phys. Chem. 1993, 97, 6859-6866 and available from:Gas Phase Core Excitation Database: http://unicorn.mcmaster.ca/corex/name-list.html
	
INFORMATION ON THE E2 (COMPARISON) CODE	

usage:  

E2.py [-h] [-ft {Athena,user_defined}] [-offset OFFSET] FILE column_energy column_intensity n_columns
											 
E2: Experimental spectra peak fitting
											 									 
positional arguments:

  FILE                  The experimental spectra file to be read in;
                        filename must include its path and it is best not to be in the directory structure for this code
  
  column_energy         The column in spectra file that holds the energy, 0 is
                        the lowest value.
						
  column_intensity      The column in spectra file that holds the intensity, 0
                        is the lowest value.
						
  n_columns             The number of columns in spectra file.

  
optional arguments:

  -h, --help            show this help message and exit
  
  -ft {Athena,user_defined}, --file_type {Athena,user_defined}
                        Select the default set of orca parameters for
                        particular chemical states.
						
  -offset OFFSET        The number of lines in the input spectra file that are
                        to be skipped before the data is read in.

EXAMPLE USAGE

1, minimal needed to run (defaults assume it is a gas):

E2.py  C:\Users\userid\Desktop\analysis\ImidazoleExperiment\N1s_Imidazole_ISEELS.txt 0 2 6

BACKGROUND

This program is the E2 part of the experimental pipeline or workflow, which is described in the README file in the top level directory for the project.
