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


exp_data_filename=r"N1s_Imidazole_ISEELS.txt" #experimental data file provided by the user
theory_data_filename=r"TDDFT_N-edge.out"
fitted_peaks_param_filename=r"fitted_peaks_param.txt"
#input and output files 
working_dir=os.getcwd()
#working_dir=r"C:\Users\dmg81179\Desktop\Code_Development\2018_June_Compare_C1\working_dir"
path_in=working_dir
path_out=working_dir
exp_data_file=path_in+r"\%s" %exp_data_filename 
edge_data_table=path_in+r"\edge_data.txt"
fitted_peaks_param=path_in+r"\%s" %fitted_peaks_param_filename
tddft_output_file=path_in+r"\%s" %theory_data_filename
theory_data_file=path_in+r"\%s.abs.dat" %theory_data_filename

###Extract Experimental data
#put experimental data into array
exp_data=np.array([])
with open (exp_data_file) as exp_data_file:
    lines_after_heading=exp_data_file.readlines()[38:]# 38 is default but it can change 
    for line in lines_after_heading:
        line=line.split()
        exp_data=np.append(exp_data,line)
    exp_data=exp_data.reshape(int(len(exp_data)/6),6)#6 is the number of columns and it can change 
exp_data_file.close()
#xdata and ydata
exp_ycol=2
exp_xdata=exp_data[:,0].astype(float) #xdata is always the first column
exp_ydata=exp_data[:,exp_ycol].astype(float) # ask user for the column where the ydata is 
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
e0_pos=np.where(dif1==max(dif1))[0][0]
e0=exp_xdata[e0_pos]
e_edge=abs(edge_data[:,4].astype(float)-e0)
e1_pos=np.where(e_edge==min(e_edge))[0][0]
element=edge_data[e1_pos][0]
#determine fitting range
fit_range=[5,5]
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
#extract Loewdin orbital population 
lines=[]
with open (tddft_output_file, "r") as tddft_output_file:
    copy=False
    for line in tddft_output_file:
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER" in line.strip():
            copy=True
        elif "MAYER POPULATION ANALYSIS" in line.strip():
            copy=False
        elif copy:
            line=line.split()
            lines.append(line)
    lines=lines[2:-3]
    tddft_output_file.close()
#This was taken from here https://codereview.stackexchange.com/questions/179530/
#split-list-of-integers-at-certain-value-efficiently
def block_split(seq,condition):
    group=[]
    for element in seq:
        if element != condition:
            group.append(element)
        elif group:
            yield group
            group = []

blocks=list(block_split(lines, []))
Loewdin_population_per=[]            
for chunk in blocks:
    states=chunk[0]   
    for k in chunk[4:-1]:
        for j in range(0,len(states)):
            Loewdin_population_per.append([int(states[j]),str(k[0]),str(k[1]),str(k[2]),float(k[j+3])])
            Loewdin_population_per=sorted(Loewdin_population_per)

#put theoretical data into array 
theory_data=np.array([])
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
fit_range=[5,5]
fit_pos_i_ls=abs(theory_xdata-(e0-12)+fit_range[0])
fit_pos_i=np.where(fit_pos_i_ls==min(fit_pos_i_ls))[0][0]
fit_pos_f_ls=abs(theory_xdata-(e0-12)-fit_range[1])
fit_pos_f=np.where(fit_pos_f_ls==min(fit_pos_f_ls))[0][0]
fit_theory_xdata=theory_xdata[fit_pos_i:fit_pos_f].astype(float)
fit_theory_ydata=norm_theory_ydata[fit_pos_i:fit_pos_f].astype(float)
#determination of number of gaussian and their initial position guesses
#smooth_ydata=scipy.signal.savgol_filter(fit_theory_ydata,11,2)
dif2=np.diff(fit_theory_ydata)/np.diff(fit_theory_xdata)
e0_pos2=np.where(dif2==max(dif2))[0][0]
posit=np.array((np.where((dif2[1:]<0)*(dif2[0:-1]>0))),dtype='int')+1
posit=np.unique(np.sort(np.append(posit,np.array((np.where((dif2[1:]>0)*(dif2[0:-1]<0))),dtype='int')+1)))
posit=posit[posit>e0_pos2]
funccenter=fit_theory_xdata[posit]
#remove extraneous peaks
b=[]
b_new=[]
for i in range(0,len(funccenter)-1):
    if funccenter[i+1]-funccenter[i]>0.5:
        b.append(i)
        b.append(i+1)
for j in b:
    if j not in b_new:
        b_new.append(j)
funccenter=funccenter[b_new]
#translate the energy scale for the theoretical data based on experimental data
trans_theory_xdata=theory_xdata+(float(first_exp_peak_center)-min(funccenter))


###plotting
fig_raw=plt.figure()
plt.plot(exp_xdata, exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(theory_xdata, theory_ydata, 'g--',label='Theoretical Spectrum')
plt.legend(loc='upper right')
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
fig_raw.savefig('Comparison of raw experimental and theoretical spectra.tif')

fig_norm=plt.figure()
plt.plot(exp_xdata, norm_exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(theory_xdata, norm_theory_ydata, 'g--',label='Theoretical Spectrum')
plt.legend(loc='upper right')
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
fig_norm.savefig('Comparison of normalized experimental and theoretical spectra.tif')

fig_trans=plt.figure()
plt.plot(exp_xdata, norm_exp_ydata, 'b',label='Experimental Spectrum')
plt.plot(trans_theory_xdata, norm_theory_ydata, 'g--',label='Theoretical Spectrum')
plt.legend(loc='upper left')
ax = plt.gca()
ax.set_xlim([fit_exp_xdata[0],fit_exp_xdata[-1]+1])
ax.set_ylim([0,max(fit_exp_ydata)+0.2])
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
fig_trans.show()
fig_trans.savefig('Comparison of normalized and translated experimental and theoretical spectra.tif.tif')


stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 