# -*- coding: utf-8 -*-
"""
Created on Tue May 22 17:16:41 2018
@author: Laila Al-Madhagi 
email: fy11lham@leeds.ac.uk 
This is a python code for peak fitting with some ideas from 
Matlab code published in: Journal of Physics: Conference Series 712 (2016) 012070  
"""
import timeit
start = timeit.default_timer()
import os
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
import scipy.signal


#exp_data_filename=r"N1s_Imidazole_ISEELS.txt" #experimental data file provided by the user
#theory_data_filename=r"TDDFT_N-edge.out"
#fitted_peaks_param_filename=r"fitted_peaks_param.txt"
Loewdin_filename=r"TDDFT_N-edge.out"
#input and output files 
working_dir=os.getcwd()
#working_dir=r"C:\Users\dmg81179\Desktop\Code_Development\2018_June_Compare_C1\working_dir"
path_in=working_dir
path_out=working_dir
#exp_data_file=path_in+r"\%s" %exp_data_filename 
#edge_data_table=path_in+r"\edge_data.txt"
#fitted_peaks_param=path_in+r"\%s" %fitted_peaks_param_filename
#tddft_output_file=path_in+r"\%s" %theory_data_filename
#theory_data_file=path_in+r"\%s.abs.dat" %theory_data_filename
Loewdin_file=path_in+r"\%s" %Loewdin_filename





lines=[]
Loewdin_pop=[]
with open (Loewdin_file, "r") as Loewdin_file:
    copy=False
    for line in Loewdin_file:
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER" in line.strip():
            copy=True
        elif "MAYER POPULATION ANALYSIS" in line.strip():
            copy=False
        elif copy:
            line=line.split()
            lines.append(line)
    lines=lines[2:-4]
    Loewdin_file.close()
"""   
states=lines[0]   
for k in lines[4:-1]:
    for j in range(0,len(states)):
        Loewdin_pop.append([states[j],k[0],k[1],k[2],k[j+3]])
        Loewdin_pop=sorted(Loewdin_pop)
"""    
    