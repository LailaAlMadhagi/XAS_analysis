XAS_analysis: SW1 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/05/17. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later

This program was developed with Python 3.5.2 with modules: subprocess and numpy(version 1.13.1). It was developed on a windows system using Anaconda.


FILES INCLUDED 

1,) README.md				This currently holds all the documentation for
							this repository

2,) SW1.py  				Python wrapper source code calculations of a gas
						
						
3,) S_args.txt			 	This text file contains all the arguements to be 
							passed from the wrapper script to the other scripts.
							
INFORMATION ON THE SW1 (SOLID WRAPPER) CODE

usage: 

LW1.py [-h] FILE


SW1: SOLID WRAPPER

positional arguments:

  FILE					Input file with arguments to be passed to other scripts 
						File name must include its path and it should be in the 
						directory structure for this code

optional arguments:

  -h, --help          	show this help message and exit
  
	

EXAMPLE USAGE

1, minimal needed to run:

GW.py C:\Users\userid\Desktop\XAS_Analysis\SW\S_args.txt


BACKGROUND

This program is the SW part of the pipeline or workflow. This is a wrapper script to run other scripts needed for solid phase calculations. The computation is in 2 steps:
1, run E2 2, run SES1 and C1 in a loop until theoretical and experimental data compare well.
