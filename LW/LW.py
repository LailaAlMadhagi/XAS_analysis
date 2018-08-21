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
import argparse
import sys
import os
import datetime
import subprocess as sp

# handle the input flags
description='LW: Wrapper script for Liquid samples'
parser = argparse.ArgumentParser(description)



parser.add_argument('in_args',
    type=argparse.FileType('r'),
    help="The arguements file to be read in.",
    default=sys.stdin, 
    metavar="FILE")

args = parser.parse_args()

path_LW, args_file = os.path.split(args.in_args.name)
LW_date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
resultsdir = r"Results_dir_"+LW_date_time

path_in=path_LW
path_out=path_LW+r'//'+resultsdir
os.makedirs(path_out)

arguments_d={'geom_directory':[],'orca_param':[],'orca_executable':[],
             'experimental_spectra':[],'experimental_energy_column_number':[],
             'experimental_intensity_column_number':[],'experimental_number_columns':[],
             'experimental_header_skip':[],'element_calculate':[],'results_dir':path_out,
             'geom_file_name':[]}

with open (path_LW+r'\%s'%args_file,'r') as args_f:
    lines=args_f.readlines()[1:]
    for line in lines:
        arguments_d[line.split('=')[0]]=line.split('=')[1].replace('\n','')
args_f.close()

#run LE2:
LE2_p=sp.Popen(['python','U:\XAS_analysis\LE2\LE2.py',
                str(arguments_d['experimental_spectra']),
                str(arguments_d['experimental_energy_column_number']),
                str(arguments_d['experimental_intensity_column_number']),
                str(arguments_d['experimental_number_columns']),
                '-offset',str(arguments_d['experimental_header_skip']), 
                '-path_out',str(arguments_d['results_dir'])],stdout=sp.PIPE, stderr=sp.PIPE)
LE2_output, LE2_err = LE2_p.communicate()
LE2_p_status=LE2_p.wait()
LE2_path_out=LE2_output.decode('uft-8').split('\r\n')[-1]