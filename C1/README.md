This is created by Laila Al-Madhagi (fy11lham@leeds.ac.uk) for her PhD project on the 2018/06/05. The code license should be CC-BY-NC. GitHub does not support this license type and it will be added later. 

This program was developed with Python 3.5.2 with modules: numpy(version1.13.1), scipy (version 0.19.1), matplotlib (version 2.0.2) and subprocess.It was developed on a windows system using Anaconda

FILES INCLUDED 
1,)README.md					This currently holds all the documentation for
								this repository
2,) C1.py   					Python source code for the comparison pipeline  
								
3,) edge_data.txt				A text file containing edge data information
								(derived from "Journal of Physics: Conference Series 712 (2016) 012070" which is taken from the X-ray Data Booklet)
4,) fitted_peaks_param.txt  	A text file containing the parameters for the peaks fitted to the experimental data. This file is used for testing   
								and development. 
5,) N1s_Imidazole_ISEELS.txt 	A text file containing experimental nitrogen K-edge data of imidazole in the gas phase. The file is used for testing and development 
								(data published in Apen et al, J. Phys. Chem. 1993, 97, 6859-6866 and available from:Gas Phase Core Excitation Database: http://unicorn.mcmaster.ca/corex/name-list.html
6,) TDDFT_N-edge.out 			TDDFT outputfile from TE1 pipeline. The file is used for testing and development
								
								
								
This program is the C1 part of the Manual pipeline, a diagram of all the planned pipes is provided in the respository.The computation is done by comparing the experimental and theoretical spectra. The theoretical spectra is translated so that the first peak is aligned with the first peak in the experimental spectrum.   