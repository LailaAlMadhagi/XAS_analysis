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
#import scipy.signal
import datetime
import argparse
import sys
import socket
from collections import OrderedDict
import shutil


# handle the input flags
description='C1: Compares experimental spectra with theorectically calculated data.'
parser = argparse.ArgumentParser(description)

parser.add_argument('in_experiment_file',
    type=argparse.FileType('r'),
    help="Experimental spectra datafile to be read in.",
    default=sys.stdin, metavar="FILE")

parser.add_argument("column_energy",  
                    type=int, 
                    help="The column in spectra file that holds the energy, 0 is the lowest value.")

parser.add_argument("column_intensity",  
                    type=int, 
                    help="The column in spectra file that holds the intensity, 0 is the lowest value.")

parser.add_argument("n_columns",  
                    type=int, 
                    help="The number of columns in spectra file.")

parser.add_argument("-offset", 
                    dest="offset",
                    type=int,
                    required=False,
                    help="The number of lines in the input spectra file that are to be skipped before the data is read in.")

parser.add_argument('in_theoretical_file',
    type=argparse.FileType('r'),
    help="Theoretically calulated spectra datafile to be read in.",
    default=sys.stdin, metavar="FILE")

parser.add_argument('in_fitted_peaks_file',
    type=argparse.FileType('r'),
    help="Peaks fitted to the experimental spectra datafile to be read in.",
    default=sys.stdin, metavar="FILE")


'''parser.add_argument('in_edge_data_file',
    type=argparse.FileType('r'),
    help="Edge data from experimental spectra datafile to be read in.",
    default=sys.stdin, metavar="FILE")'''


args = parser.parse_args()

skip = int(args.offset)

file_exp_read=args.in_experiment_file.name

path, file_exp = os.path.split(file_exp_read)
                   
index_of_dot = file_exp.index(".") 
file_exp_without_extension = file_exp[:index_of_dot]

date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
resultsdir = r"C1_"+file_exp_without_extension+r"_"+date_time


path_in=path
path_out=path+r"/"+resultsdir
os.makedirs(path_out)

log_file_name = path_out+r"\\log.txt"
log_file=open(log_file_name, "w") 

log_file.write(description+"\n\n")
host=socket.gethostbyaddr(socket.gethostname())[0]
log_file.write(r"This program ran at "+date_time+r" on the "+host+r" host system.")
log_file.write("\n\n")

html_line1 =r"This program ran at "+date_time+r" on the "+host+r" host system"




print("\n~ Experimental spectra file details: {}".format(args.in_experiment_file))
log_file.write("\n\n~ Experimental spectra file details: {}".format(args.in_experiment_file))


print("\n~ Theoretical spectra file details: {}".format(args.in_theoretical_file))
log_file.write("\n\n~ Theoretical spectra file details: {}".format(args.in_theoretical_file))


print("\n~ Peaks fitted to the theoretical spectra file details: {}".format(args.in_fitted_peaks_file))
log_file.write("\n\n~ Peaks fitted to the theoretical spectra file details: {}".format(args.in_fitted_peaks_file))




#input and output files 
working_dir=os.getcwd()


#exp_data_file=path_in+r"\%s" %exp_data_filename 
exp_data_file=args.in_experiment_file.name
#edge_data_table=path_in+r"\edge_data.txt"
#edge_data_table=args.in_edge_data_file.name
edge_data_table=working_dir+r"\..\edge_data.txt"
#fitted_peaks_param=path_in+r"\%s" %fitted_peaks_param_filename
fitted_peaks_param=args.in_fitted_peaks_file.name
#theory_data_file=path_in+r"\%s.abs.dat" %theory_data_filename
file_theory_read=args.in_theoretical_file.name
path_theory, file_theory = os.path.split(file_theory_read)
theory_data_file=path_out+"/%s.abs.dat" %file_theory
theory_data_filename=path_out+r"/%s" %file_theory
tddft_output_file=r"%s" %theory_data_filename

html_infile_name=working_dir+r"\template.html"
html_outfile_name=path_out+r"\report.html"



        
        



print("theory_data_file: ",theory_data_file)


###move theory data file to results directory
path, theoretical_file_name=os.path.split(args.in_theoretical_file.name)
shutil.copy(theoretical_file_name,path_out)

###Extract Experimental data
#put experimental data into array
exp_data=np.array([])
with open (exp_data_file) as exp_data_file:
    lines_after_heading=exp_data_file.readlines()[skip:] 
    for line in lines_after_heading:
        line=line.split()
        exp_data=np.append(exp_data,line)
    exp_data=exp_data.reshape(int(len(exp_data)/args.n_columns),args.n_columns)
exp_data_file.close()
#xdata and ydata
exp_ycol=args.column_intensity
#exp_xdata=exp_data[:,args.column_energy].astype(float) #xdata is always the first column
exp_xdata_row=exp_data[:,args.column_energy]
exp_xdata=np.array([])
for x in exp_xdata_row:
    exp_xdata=np.append(exp_xdata,float(x.replace(",","")))
#exp_ydata=exp_data[:,exp_ycol].astype(float) # ask user for the column where the ydata is 
exp_ydata_row=exp_data[:,exp_ycol]
exp_ydata=np.array([])
for y in exp_ydata_row:
    exp_ydata=np.append(exp_ydata,float(y.replace(",","")))
norm_exp_ydata=(exp_ydata-min(exp_ydata))/(max(exp_ydata)-min(exp_ydata))

#put edge data into array
edge_data=np.array([])
with open(edge_data_table,"r") as edge_data_table:
    next (edge_data_table)
    next (edge_data_table)
    for line in edge_data_table:
        line=line.split()
        edge_data=np.append(edge_data,line)
    edge_data=edge_data.reshape(int(len(edge_data)/5),5)
    edge_data_table.close()
##determine fitting range to be used when plotting 
#determine e0 requried for fitting range
dif1=np.diff(exp_ydata)/np.diff(exp_xdata)
#e0_pos=np.where(dif1==max(dif1))[0][0]
e0_pos=np.where(dif1==np.nanmax(dif1[dif1!=np.inf]))[0][0]
e0=exp_xdata[e0_pos]
e_edge=abs(edge_data[:,4].astype(float)-e0)
e1_pos=np.where(e_edge==min(e_edge))[0][0]
element=edge_data[e1_pos][0]
#determine fitting range
fit_range=[5,8]
fit_pos_i_ls=abs(exp_xdata-e0+fit_range[0])
fit_pos_i=np.where(fit_pos_i_ls==min(fit_pos_i_ls))[0][0]
fit_pos_f_ls=abs(exp_xdata-e0-fit_range[1])
fit_pos_f=np.where(fit_pos_f_ls==min(fit_pos_f_ls))[0][0]
fit_exp_xdata=exp_xdata[fit_pos_i:fit_pos_f].astype(float)
fit_exp_ydata=norm_exp_ydata[fit_pos_i:fit_pos_f].astype(float)

#read fitted peaks parameter 
peak_param=[]
with open(fitted_peaks_param) as fit_param:
    for line in fit_param:
       line=line.split()
       peak_param.append(line)
    fit_param.close()

for element in peak_param:
    if "center" in element:
        first_exp_peak_center=min(element)
###

###Extract theoretical data
#run tddft calc
p=sp.Popen(['orca_mapspc',tddft_output_file,'ABS','-eV','-x0%s'%str(int(min(fit_exp_xdata)-15)),'-x1%s'%str(int(max(fit_exp_xdata))),'-n500','-w0.6'])
p_status=p.wait()

#read orbital energy information from ouput file
#might not this information now
""" 
orbital_energies_array=np.array([])
with open(theory_data_file, 'r+') as tddft_output_file:
    copy=False
    for line in tddft_output_file:
        if line.strip()=="ORBITAL ENERGIES":
            copy=True
        elif line.strip()=="MOLECULAR ORBITALS":
            copy=False
        elif copy:
            line=line.split()
            if len(line)==4:
                orbital_energies_array=np.append(orbital_energies_array,line)
    orbital_energies_array=orbital_energies_array.reshape(int(len(orbital_energies_array)/4),4)
    orbital_energies_array=np.delete(orbital_energies_array, np.s_[0:int(len(orbital_energies_array)/2)],0)             
    unocc_orbitals=[]
    for x in range (0, int(len(orbital_energies_array)-1)):
        if orbital_energies_array[x][1]==np.str("0.0000"):
            unocc_orbitals.append(x)
        x+=1
    orbital_energies_array=np.delete(orbital_energies_array, np.s_[unocc_orbitals[0]:],0)
    tddft_output_file.close()
"""
#extract Loewdin orbital population analysis information and excited states information
Loewdin_lines=[]
states_lines=[]
with open (tddft_output_file, "r") as tddft_output_file:
    copy=False
    for line in tddft_output_file:
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER" in line.strip():
            copy=True
        elif "MAYER POPULATION ANALYSIS" in line.strip():
            copy=False
        elif copy:
            line=line.split()
            Loewdin_lines.append(line)
    tddft_output_file.seek(0)
    copy=False
    for line in tddft_output_file:
        if "EXCITED STATES" in line.strip():
            copy=True
        elif "EXCITATION SPECTRA" in line.strip():
            copy=False
        elif copy:
            line=line.split()
            states_lines.append(line)
    
    tddft_output_file.close()
#Function to split read lines into blocks
#The function was adapted from here https://codereview.stackexchange.com/questions/179530/
#split-list-of-integers-at-certain-value-efficiently
def block_split(seq,condition):
    group=[]
    for element in seq:
        if condition not in element and element!=[]:
            group.append(element)
        elif group:
            yield group
            group = []
def states_split(seq,condition):
    group=[]
    for element in seq:
        if condition in element and element!=[]:
            group.append(element)
        elif group:
            yield group
            group = []

Loewdin_lines=Loewdin_lines[2:-3]
Loewdin_blocks=list(block_split(Loewdin_lines, []))
Loewdin_population_per=[]            
for chunk in Loewdin_blocks:
    states=chunk[0]   
    for k in chunk[4:-1]:
        for j in range(0,len(states)):
            Loewdin_population_per.append([int(states[j]),str(k[0]),str(k[1]),str(k[2]),float(k[j+3])])
            Loewdin_population_per=sorted(Loewdin_population_per)
            
states_lines=states_lines[4:-1]
states=list(states_split(states_lines,'STATE'))
states_blocks=list(block_split(states_lines, 'STATE'))
excited_states_ls=[]
for state,chunk in zip(states,states_blocks):
    for i in chunk:
        excited_states_ls.append([int(state[0][1].replace(":","")),float(state[0][5]),i[0],i[2],float(i[4])])                        

#put theoretical data into array 
theory_data=np.array([])
#with open (theory_data_file) as theory_data_file:
with open (theory_data_file) as theory_data_file:
    for line in theory_data_file:
        line=line.split()
        theory_data=np.append(theory_data,line)
    theory_data=theory_data.reshape(int(len(theory_data)/5),5)#6 is the number of columns and it can change 
theory_data_file.close()
#xdata and ydata
theory_ycol=1
theory_xdata=theory_data[:,0].astype(float) 
theory_ydata=theory_data[:,theory_ycol].astype(float) 
norm_theory_ydata=(theory_ydata-min(theory_ydata))/(max(theory_ydata)-min(theory_ydata))

##this is done to determine the center of the first peak in the theoretical spectrum
#determine fitting range
fit_range=[5,8]
fit_pos_i_ls=abs(theory_xdata-(e0-12)+fit_range[0])
fit_pos_i=np.where(fit_pos_i_ls==min(fit_pos_i_ls))[0][0]
fit_pos_f_ls=abs(theory_xdata-(e0-12)-fit_range[1])
fit_pos_f=np.where(fit_pos_f_ls==min(fit_pos_f_ls))[0][0]
fit_theory_xdata=theory_xdata[fit_pos_i:fit_pos_f].astype(float)
fit_theory_ydata=norm_theory_ydata[fit_pos_i:fit_pos_f].astype(float)
#determination of number of gaussian and their initial position guesses
#smooth_ydata=scipy.signal.savgol_filter(fit_theory_ydata,11,2)
dif2=np.diff(fit_theory_ydata)/np.diff(fit_theory_xdata)
e0_pos2=next((i for i, x in enumerate(dif2) if x), None)
posit=np.array((np.where((dif2[1:]<0)*(dif2[0:-1]>0))),dtype='int')+1
posit=np.unique(np.sort(np.append(posit,np.array((np.where((dif2[1:]>0)*(dif2[0:-1]<0))),dtype='int')+1)))
posit=posit[posit>e0_pos2]
funccenter=fit_theory_xdata[posit]
#remove extraneous peaks
b=[]
b_new=[]
for i in range(0,len(funccenter)-1):
    if funccenter[i+1]-funccenter[i]>0.2:
        b.append(i)
        b.append(i+1)
for j in b:
    if j not in b_new:
        b_new.append(j)
funccenter=funccenter[b_new]


###Peak assignement
peak_assignment_ls=[]
for j in excited_states_ls:
    for k in Loewdin_population_per:
        if str(k[0]) in j[2]:
            if k[3]=='s':
                if k[4]>=5:
                    j.append(k[1]+k[2])                     
for i in funccenter:
    for j in excited_states_ls:
        if j[1]-0.5 <= i <= j[1]+0.5:
            for k in Loewdin_population_per:
                if str(k[0]) in j[3]:
                    if len(j)>5:
                        if k[1]+k[2] in j[5]:
                            if k[3]=='pz':
                                if k[4]>=5:
                                    peak_assignment_ls.append([float("%.2f"%i),j[1],j[2],j[5],j[3],j[4]])
                    else:
                        pass

peak_assignment_ls_f=[]
if len(peak_assignment_ls)>3:
    for index, element in enumerate(peak_assignment_ls):
        if element[1]==peak_assignment_ls[index-2][1]:
                pass
        else:
            peak_assignment_ls_f.append(element)
else:
    peak_assignment_ls_f=peak_assignment_ls
        
peak_assignment_d=OrderedDict()
for element in peak_assignment_ls_f:
    if element[0] not in peak_assignment_d:
        peak_assignment_d[element[0]]=[element[1:]]
    else:
        peak_assignment_d[element[0]].append(element[1:])

###
#translate the energy scale for the theoretical data based on experimental data
if peak_assignment_ls_f!=[]:
    transform=float(first_exp_peak_center)-peak_assignment_ls_f[0][0]
else:
    transform=float(first_exp_peak_center)-min(funccenter)
trans_theory_xdata=theory_xdata+transform
html_transform="%.3f" % transform

###plotting
fig_raw=plt.figure(num=None, figsize=(10, 8), dpi=600, facecolor='w', edgecolor='k')
plt.plot(exp_xdata, exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(theory_xdata, theory_ydata, 'g--',label='Theoretical Spectrum')
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig_raw.savefig(path_out+r'\\ComparisonOfRawExperimentalAndTheoreticalSpectra.png',bbox_inches='tight')

fig_norm=plt.figure(num=None, figsize=(10, 8), dpi=600, facecolor='w', edgecolor='k')
plt.plot(exp_xdata, norm_exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(theory_xdata, norm_theory_ydata, 'g--',label='Theoretical Spectrum')
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig_norm.savefig(path_out+r'\\ComparisonOfNormalizedExperimentalAndTheoreticalSpectra.png',bbox_inches='tight')

fig_trans=plt.figure(num=None, figsize=(10, 8), dpi=600, facecolor='w', edgecolor='k')
plt.plot(exp_xdata, norm_exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(trans_theory_xdata, norm_theory_ydata, 'g--',label='Theoretical Spectrum')
ax = plt.gca()
ax.set_xlim([fit_exp_xdata[0],fit_exp_xdata[-1]+1])
ax.set_ylim([0,max(fit_exp_ydata)+0.2])
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig_trans.show()
fig_trans.savefig(path_out+r'\\ComparisonOfNormalizedAndTranslatedExperimentalAndTheoreticalSpectra.png',bbox_inches='tight')

fig_trans=plt.figure(num=None, figsize=(10, 8), dpi=600, facecolor='w', edgecolor='k')
plt.plot(exp_xdata, norm_exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(trans_theory_xdata, norm_theory_ydata, 'g--',label='Theoretical Spectrum')
ax = plt.gca()
ax.set_xlim([fit_exp_xdata[0],fit_exp_xdata[-1]+1])
ax.set_ylim([0,max(fit_exp_ydata)+0.2])
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
feature=1
for element in peak_assignment_d:
    plt.axvspan(element+transform,element+transform, facecolor='g', alpha=1)
    peak_assign_string=str(feature)
    #for item in peak_assignment_d[element]:
        #peak_assign_string+='%s (%s > %s) (%s)\n' %(item[2],item[1],item[3],round(item[4],3))
    plt.annotate(peak_assign_string, 
             xy=(element+transform, float(norm_theory_ydata[np.where(np.around(trans_theory_xdata,6)==round(element+transform,6))])),
             arrowprops=dict(arrowstyle="->", connectionstyle="angle3",lw=1))
    feature+=1
fig_trans.show()
fig_trans.savefig(path_out+r'\\ComparisonOfNormalizedAndTranslatedExperimentalAndTheoreticalSpectraWithPeakAssignment.png',bbox_inches='tight')

stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 

log_file.write("\n\nEND:\nRunning time is: "+ str(round(running_time,3)) + " minutes")

log_file.close()


# The resolution of the images depends on the graphics hardware and settings.
# The html template scales the images so the image resolution is not so important
# Image scaling is not part of html 5 format and at the time of writing html 5
# was not working in all browsers so the template is not in html 5.
# These are also described in the README file.

feature=1

with open(html_infile_name, "r") as html_in, open(html_outfile_name, "w") as html_out:
    n=0
    for line in html_in:
        if '***' in line and n==2:
            html_out.write(line.replace("***",html_transform))            
            n+=1
        elif '***' in line and n==1:
            html_out.write(line.replace("***",html_line1))            
            n+=1
        elif '***' in line and n==0:
            html_out.write(line.replace("***",description))
            n+=1
        elif '+*+*+*' in line:
            #print("before: ",line)
            line=""
            #html_out.write(line.replace("+*+*+*","\n\n"))
            print("after: ",line)
            #s0=1
            for element in peak_assignment_d:
                #print("1: ",element)                
                for item in peak_assignment_d[element]:
                    peak_assign_string=""
                    peak_assign_string+='<tr> <td>%s</td> <td>%s</td> <td>%s</td> <td>%s</td> <td>%s</td> </tr>\n' %(feature,item[2],item[1],item[3],round(item[4],3))
                    html_out.write(peak_assign_string)
                    #print("2: ",peak_assign_string)
                feature+=1
        elif '+++' in line:
            line=""
            s0=0
            for index, element in enumerate(peak_assignment_ls_f):
                peak_number=0
                if element[1]==peak_assignment_ls_f[index-1][1]:
                    peak_number=str(s0)
                else:
                    s0+=1
                    peak_number=str(s0)
                s1="%.3f" % element[1]
                new_line=""
                new_line='<tr> <td> %d </td> <td> %s </td> <td> %s </td> <td> %s </td> <td> %s(%f) </td> </tr> \n' %(s0,s1,element[2],element[3],element[4],round(element[5],3))
                html_out.write(new_line)
        else:
            html_out.write(line.strip())
        

            
