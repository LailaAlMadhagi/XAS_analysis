XAS_analysis: GC1 code

GENERAL INFORMATION ON THIS PROJECT
This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/06/05. Joanna Leng (j.leng@leeds.ac.uk) from the University of Leeds contributed to the development. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later. 

This program was developed with Python 3.5.2 with modules: numpy(version1.13.1), scipy (version 0.19.1), matplotlib (version 2.0.2) and subprocess.It was developed on a windows system using Anaconda

FILES INCLUDED 

1,) README.md					This currently holds all the documentation for
								this repository
								
2,) C1.py   					Python source code for the comparison pipeline  
								
3,) fitted_peaks_param.txt  	A text file containing the parameters for
								the peaks fitted to the experimental data. This file is used for testing   
								and development. 
								
5,) N1s_Imidazole_ISEELS.txt 	A text file containing experimental nitrogen K-edge data of imidazole in the gas phase. The file is used for testing and development 
								(data published in Apen et al, J. Phys. Chem. 1993, 97, 6859-6866 and available from:Gas Phase Core Excitation Database: http://unicorn.mcmaster.ca/corex/name-list.html
								
6,) TDDFT_N-edge.out 			TDDFT outputfile from TE1 pipeline. The
								file is used for testing and development.
								
								
								
INFORMATION ON THE GC1 (GAS COMPARISON) CODE

usage:   

[-h] [-offset] [-orca] [-path_out] FILE column_energy column_intensity n_columns FILE FILE

GC1: Compares experimental spectra with theorectically calculated data for a gas.
	   
	   
positional arguments:

  FILE        			Experimental spectra datafile to be read in.
  
  column_energy         The column in spectra file that holds the energy, 0 is
                        the lowest value.
						
  column_intensity      The column in spectra file that holds the intensity, 0
                        is the lowest value.
						
  n_columns             The number of columns in spectra file.
  
  FILE        			Theoretically calulated spectra datafile to be read in.
  
  FILE        			Peaks fitted to the experimental spectra datafile to be read in.

optional arguments:

  -h, --help  			show this help message and exit
  
  -offset OFFSET        The number of lines in the input spectra file that are
                        to be skipped before the data is read in.
						
  -orca FILE          	path to the orca executable; C:\Orca\orca is the
						default path-orca	
  
  -path_out				The directory where the results directory will be written
  
EXAMPLE USAGE

C1.py C:\Users\userid\Desktop\data\Imidazole\N1s_Imidazole_ISEELS.txt 0 2 6 -offset 38 C:\Users\userid\Desktop\data\Imidazole\TDDFT_N-edge.out C:\Users\userid\Desktop\data\Imidazole\N1s_Imidazole_ISEELS.txt_fitted_peaks_param.txt

THE REPORT

Images in this script are created and are referenced in an html report. The images exist separately to this report and can be used in other documents. The resolution of these images is determined by the graphics settings of the system that creates them. These are scaled in the html to 100% so that high resolution images will be visible on the screen (it will also increase the size of low resolution images). At the time of writing the html 5 format was new and was not fully implemented in all browsers. Also the html 5 standard had new ways to handle images that did not include scaling so the html template is not compatible with html 5.

  
BACKGROUND  

This program is the GC1 part of the Manual pipeline, which is explained the README file in the top directory of this project.The computation is done by comparing the experimental and theoretical spectra. The theoretical spectra is translated so that the first peak is aligned with the first peak in the experimental spectrum.   
