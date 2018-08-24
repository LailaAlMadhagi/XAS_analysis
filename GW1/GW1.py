# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 00:46:46 2018
@author: Laila Al-Madhagi 
email: fy11lham@leeds.ac.uk 
"""
import timeit
start = timeit.default_timer()
import argparse
import sys
import os
import datetime
import subprocess as sp
import numpy as np


# handle the input flags
description='GW: Wrapper script for gas samples'
parser = argparse.ArgumentParser(description)



parser.add_argument('in_args',
    type=argparse.FileType('r'),
    help="The arguements file to be read in.",
    default=sys.stdin, 
    metavar="FILE")

args = parser.parse_args()

path_args, args_file = os.path.split(args.in_args.name)
working_dir=os.getcwd()
parent_dir=path_args+r'//..'
GW_date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
resultsdir = r"Results_dir_"+GW_date_time

path_in=working_dir
path_out=path_in+r'//'+resultsdir
os.makedirs(path_out)

# set log file
log_file_name = path_out+r"//GW_log.txt"
log_file=open(log_file_name, "w") 
log_file.write(description+"\n\n")

arguments_d={'geom_file_name':[],'orca_param':[],'orca_executable':[],
             'experimental_spectra':[],'experimental_energy_column_number':[],
             'experimental_intensity_column_number':[],'experimental_number_columns':[],
             'experimental_header_skip':[],'element_calculate':[],'results_dir':path_out,
             'tddft_out_file':[],'fitted_peaks_params':[]}

with open (args.in_args.name,'r') as args_f:
    lines=args_f.readlines()[1:]
    for line in lines:
        arguments_d[line.split('=')[0]]=line.split('=')[1].replace('\n','')
args_f.close()

#run E2:
E2_p=sp.Popen(['python',str(parent_dir+r'//E2//E2.py'),
                str(arguments_d['experimental_spectra']),
                str(arguments_d['experimental_energy_column_number']),
                str(arguments_d['experimental_intensity_column_number']),
                str(arguments_d['experimental_number_columns']),
                '-offset',str(arguments_d['experimental_header_skip']), 
                '-path_out',str(arguments_d['results_dir'])],stdout=sp.PIPE, stderr=sp.PIPE)
E2_output, E2_err = E2_p.communicate()
E2_p_status=E2_p.wait()
E2_path_out=((E2_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]
log_file.write('E2 output is: '+E2_output.decode('utf-8'))
log_file.write("\n\n")
log_file.write('E2 err is: '+E2_err.decode('utf-8'))
log_file.write("\n\n")
log_file.flush()

R_sqr=0
theory_xdata_all=np.array([])
theory_xdata_avg=np.array([])
theory_ydata_all=np.array([])
theory_ydata_avg=np.array([])



#run GTE1
GTE1_p=sp.Popen(['python',str(parent_dir+r'//GTE1//GTE1.py'),
                str(arguments_d['geom_file_name']),
                '-orca',str(arguments_d['orca_executable']),
                '-path_out',str(arguments_d['results_dir']),
                '-element',str(arguments_d['element_calculate'])],stdout=sp.PIPE, stderr=sp.PIPE)
GTE1_output, GTE1_err = GTE1_p.communicate()
GTE1_p_status=GTE1_p.wait()
GTE1_path_out=((GTE1_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]
log_file.write('GTE1 output is: '+GTE1_output.decode('utf-8'))
log_file.write("\n\n")
log_file.write('GTE1 err is: '+GTE1_err.decode('utf-8'))
log_file.write("\n\n")
log_file.flush()
#add tddft_out file to arguments dictionary
for file in os.listdir(GTE1_path_out):
    if file.endswith('.out') and arguments_d['element_calculate'] in file:
        arguments_d['tddft_out_file']=os.path.join(GTE1_path_out,file)
#add peak_params file to arguments dirctionary
for file in os.listdir(E2_path_out):
    if file.endswith('peak_params.txt'): 
        arguments_d['fitted_peaks_params']=os.path.join(E2_path_out,file)


#run C1
C1_p=sp.Popen(['python',str(parent_dir+r'//C1//C1.py'),
                str(arguments_d['experimental_spectra']),
                str(arguments_d['experimental_energy_column_number']),
                str(arguments_d['experimental_intensity_column_number']),
                str(arguments_d['experimental_number_columns']),
                '-offset',str(arguments_d['experimental_header_skip']),
                '-orca',str(arguments_d['orca_executable']),
                str(arguments_d['tddft_out_file']),
                str(arguments_d['fitted_peaks_params']),
                '-path_out',str(arguments_d['results_dir'])],stdout=sp.PIPE, stderr=sp.PIPE)
C1_output, C1_err = C1_p.communicate()
C1_p_status=C1_p.wait()
C1_path_out=((C1_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]
log_file.write('C1 output is: '+C1_output.decode('utf-8'))
log_file.write("\n\n")
log_file.write('C1 err is: '+C1_err.decode('utf-8'))
log_file.write("\n\n")
log_file.close()            
            
            