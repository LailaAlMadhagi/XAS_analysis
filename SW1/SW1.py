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
import socket
import subprocess as sp
import numpy as np
import pandas as pd


# handle the input flags
description='SW1: Wrapper script for Solid samples'
parser = argparse.ArgumentParser(description)



parser.add_argument('in_args',
    type=argparse.FileType('r'),
    help="The arguements file to be read in.",
    default=sys.stdin, 
    metavar="FILE")

args = parser.parse_args()
#dictionary containing arguments 
arguments_d={'geom_directory':[],'orca_param':[],'orca_executable':[],
             'experimental_spectra':[],'experimental_energy_column_number':[],
             'experimental_intensity_column_number':[],'experimental_number_columns':[],
             'experimental_header_skip':[],'element_calculate':[],'results_dir':[],
             'geom_file_name':[],'tddft_out_file':[],'fitted_peaks_params':[],
             'hydrogen_opt':[],'number_of_processors':[]}

with open (args.in_args.name,'r') as args_f:
    lines=args_f.readlines()[1:]
    for line in lines:
        arguments_d[line.split('=')[0]]=line.split('=')[1].replace('\n','')
args_f.close()

path_args, args_file = os.path.split(args.in_args.name)
working_dir=os.getcwd()
script_path=os.path.dirname(os.path.abspath(__file__))    
SW1_date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
resultsdir = r"Results_dir_"+arguments_d['geom_directory'].split('/')[-1]+'_element_'+arguments_d['element_calculate']+'_'+SW1_date_time

path_in=working_dir
path_out=path_in+r'//'+resultsdir
os.makedirs(path_out)
arguments_d['results_dir']=path_out

# set log file
log_file_name = path_out+r"//SW1_log.txt"
log_file=open(log_file_name, "w") 
print(description)
log_file.write(description+"\n\n")
try:
    host=socket.gethostbyaddr(socket.gethostname())[0]
except socket.herror:
    host=''
log_file.write(r"This program ran at "+SW1_date_time+r" on the "+host+r" host system.")
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
try:
    # doesn't even have to be reachable
    s.connect(('10.255.255.255', 1))
    IP = s.getsockname()[0]
except:
    IP = '127.0.0.1'
finally:
    s.close()
    
print("IP: ",IP)
log_file.write("\nSystem's IP address is: "+IP)
log_file.write("\n\n")
log_file.flush()


#run E2:
E2_p=sp.Popen(['python',str(script_path+r'//..//E2//E2.py'),
                str(arguments_d['experimental_spectra']),
                str(arguments_d['experimental_energy_column_number']),
                str(arguments_d['experimental_intensity_column_number']),
                str(arguments_d['experimental_number_columns']),
                '-offset',str(arguments_d['experimental_header_skip']), 
                '-path_out',str(arguments_d['results_dir'])],stdout=sp.PIPE, stderr=sp.PIPE)
E2_output, E2_err = E2_p.communicate()
E2_p_status=E2_p.wait()
log_file.write('E2 output is: '+E2_output.decode('utf-8'))
log_file.write("\n\n")
log_file.write('E2 err is: '+E2_err.decode('utf-8'))
log_file.write("\n\n")
log_file.flush()
if "E2 Process Ended Successfully" not in E2_output.decode('utf-8'):
    sys.exit("Error calling running E2 script. E2 error is: "+E2_err.decode('utf-8'))
E2_path_out=((E2_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]


R_sqr=0
theory_xdata_all=np.array([])
theory_xdata_avg=np.array([])
theory_ydata_all=np.array([])
theory_ydata_avg=np.array([])

geom_filenames=[]
for file in os.listdir(arguments_d['geom_directory']):
    if not file.startswith('.') and file.endswith('.xyz'):
        geom_filenames.append(file)

loop=int(1)
loop_break=int(3)
while loop <loop_break:
    if R_sqr<0.8:
        if loop == 1:
            arguments_d['hydrogen_opt']=None
            arguments_d['geom_file_name']=geom_filenames[0]
        else:
            arguments_d['hydrogen_opt']=str(loop)  
            arguments_d['geom_file_name']=geom_filenames[0]#NEED TO change so that it is the opt geom
        #run SES1
        SES1_p=sp.Popen(['python',str(script_path+r'//..//SES1//SES1.py'),
                        str(arguments_d['geom_directory']),
                        '-geom_file_name',str(arguments_d['geom_file_name']),
                        '-orca',str(arguments_d['orca_executable']),
                        '-path_out',str(arguments_d['results_dir']),
                        '-element',str(arguments_d['element_calculate']),
                        '-h_opt',str(arguments_d['hydrogen_opt']),
                        '-pal',str(arguments_d['number_of_processors'])],stdout=sp.PIPE, stderr=sp.PIPE)
        SES1_output, SES1_err = SES1_p.communicate()
        SES1_p_status=SES1_p.wait()
        log_file.write('SES1 output is: '+SES1_output.decode('utf-8'))
        log_file.write("\n\n")
        log_file.write('SES1 err is: '+SES1_err.decode('utf-8'))
        log_file.write("\n\n")
        log_file.flush()
        if "SES1 Process Ended Successfully" not in SES1_output.decode('utf-8'):
            sys.exit("Error calling running SES1 script. SES1 error is: "+SES1_err.decode('utf-8'))
        SES1_path_out=((SES1_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]
        
        #add tddft_out file to arguments dictionary
        for file in os.listdir(SES1_path_out):
            if file.endswith('.out') and arguments_d['element_calculate'] in file:
                arguments_d['tddft_out_file']=os.path.join(SES1_path_out,file)
        #add peak_params file to arguments dirctionary
        for file in os.listdir(E2_path_out):
            if file.endswith('peak_params.txt'): 
                arguments_d['fitted_peaks_params']=os.path.join(E2_path_out,file)
        #run C1
        C1_p=sp.Popen(['python',str(script_path+r'//..//C1//C1.py'),
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
        log_file.write('C1 output is: '+C1_output.decode('utf-8'))
        log_file.write("\n\n")
        log_file.write('C1 err is: '+C1_err.decode('utf-8'))
        log_file.write("\n\n")
        log_file.flush()
        if "C1 Process Ended Successfully" not in C1_output.decode('utf-8'):
            sys.exit("Error calling running C1 script. C1 error is: "+C1_err.decode('utf-8'))
        C1_path_out=((C1_output.decode('utf-8').split('\n')[-2]).replace('\r','').replace('\n','')).split(': ')[1]

        ##comparison
        #extract normalized translated theoretical data
        for file in os.listdir(C1_path_out):
            if file.endswith('NormTranslatedTheoryData.txt'):
                Norm_trans_theory_file=os.path.join(C1_path_out,file)
        theory_xdata_s=[]
        theory_ydata_s=[]
        with open (Norm_trans_theory_file,'r') as theory_file:
            lines=theory_file.readlines()
            for line in lines:
                theory_xdata_s.append(float(line.split('\t')[0]))
                theory_ydata_s.append(float((line.split('\t')[1]).replace('\n','')))
        theory_file.close()
        theory_xdata_all=np.append(theory_xdata_all,theory_xdata_s)
        theory_ydata_all=np.append(theory_ydata_all,theory_ydata_s)
        df_theory=pd.DataFrame({'energy_theory':theory_xdata_s,'intensity_theory':theory_ydata_s}) 
        #extract experimental data needed for comparison
        exp_data=np.array([])
        with open (arguments_d['experimental_spectra'],'r') as exp_file:
            lines=exp_file.readlines()[int(arguments_d['experimental_header_skip']):]
            for line in lines:
                exp_data=np.append(exp_data,line.split())
            exp_data=exp_data.reshape(int(len(exp_data)/int(arguments_d['experimental_number_columns'])),int(arguments_d['experimental_number_columns']))
        exp_file.close()
        exp_xdata_s=exp_data[:,int(arguments_d['experimental_energy_column_number'])]
        exp_xdata=np.array([])
        for x in exp_xdata_s:
            exp_xdata=np.append(exp_xdata,float(x.replace(",","")))
        exp_ydata_s=exp_data[:,int(arguments_d['experimental_intensity_column_number'])]
        exp_ydata=np.array([])
        for y in exp_ydata_s:
            exp_ydata=np.append(exp_ydata,float(y.replace(",","")))
        norm_exp_ydata=(exp_ydata-min(exp_ydata))/(max(exp_ydata)-min(exp_ydata))
        df_exp=pd.DataFrame({'energy_exp':exp_xdata,'intensity_exp':norm_exp_ydata})
        df_C=(pd.concat([df_exp.rename(columns={'energy_exp':'energy'}),
                         df_theory.rename(columns={'energy_theory':'energy'})]).sort_values('energy')).fillna(0)
        ymean=np.sum((df_C['intensity_exp'].astype(float))/len(df_C['intensity_exp'].astype(float)))
        ssreg=np.sum((df_C['intensity_theory'].astype(float)-ymean)**2)
        sstot=np.sum((df_C['intensity_exp'].astype(float)-ymean)**2)
        R_sqr=ssreg/sstot
        loop+=1
    else:
        print("experimental and theoretical spectra are in good agreement and R-squared value is: %s"%R_sqr)
        break
if loop==loop_break:
    print("{0} hydrogen position optimization runs have been attempted. The R-squared value is: {1}".format(loop-2,'%.3f'%R_sqr))
    log_file.write("{0} hydrogen position optimization runs have been attempted. The R-squared value is: {1}".format(loop-2,'%.3f'%R_sqr))

log_file.close()            
            
            