# -*- coding: utf-8 -*-
"""
Updated on Thu May 17 10:17:57 2018 by Laila
Written by: Laila Al-Madhagi on May 11 2018
fy11lham@leeds.ac.uk
Modified by Joanna Leng on May 16 2018
j.leng@leeds.ac.uk
"""

import timeit
start = timeit.default_timer()
import subprocess as sp
import mmap
import numpy as np
import os
import datetime

print("start")

ORCA=r"C:\Orca\orca.exe"







working_dir=os.getcwd()

file_geom=""


for file_name in os.listdir(working_dir):
    if file_name.endswith('.xyz'): 
        file_geom=file_name
        print(file_geom)
               
index_of_dot = file_geom.index(".") 
file_geom_without_extension = file_geom[:index_of_dot]


resultsdir = r"TE1_"+file_geom_without_extension+r"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

print(resultsdir)
#os.makedirs(mydir)

path_in=working_dir
path_out=working_dir+r"\\"+resultsdir
os.makedirs(path_out)

# The is all the files we need to create
geom_file=path_in+r"\geom.xyz"
opt_input_file=path_out+r"\opt.inp"
new_opt_input_file=path_out+r"\new_opt.inp"
opt_output_file=path_out+r"\opt.out"
new_opt_output_file=path_out+r"\new_Opt.out"
opt_geom_file=path_out+r"\opt.xyz"
opt_2_file=path_out+r"\opt_2.xyz"
freq_input_file=path_out+r"\freq.inp"
freq_output_file=path_out+r"\freq.out"
tddft_input_file=path_out+r"\TDDFT.inp"
tddft_output_file=path_out+r"\TDDFT.out"
opt_error_file=path_out+r"\opt_error.txt"
freq_error_file=path_out+r"\freq_error.txt"
tddft_error_file=path_out+r"\TDDFT_error.txt"

orbital_energies_array=np.array([])
orbital_window_array=[]
edge_data_array=np.array([['C','6','K','0','284.2'],['N','7','K','0','409.9']])

#generate input file for Opt calculation
opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
print_array=np.array(["\n!NormalPrint","%output","Print[P_Basis] 2","Print[P_MOS] 1","end"])
geom_array=np.array(["*xyzfile","0","1",str(geom_file)])

with open(opt_input_file, "w") as opt_file:
    for item in opt_keywords_array:
        opt_file.writelines(["%s " %item])
    for item in print_array:
        opt_file.writelines(["%s \n" %item])
    for item in geom_array:
        opt_file.writelines(["%s " %item])
opt_file.close()

'''
# run Opt calculation
opt_out=open(opt_output_file, "w") 
opt_err=open(opt_error_file,"w") 
p=sp.Popen(['ORCA',opt_input_file], stdout=opt_out, stderr=opt_err)
p_status=p.wait()
opt_out.close()
opt_err.close()

# check Opt.out file
finding=-1
while (finding==-1):
    with open(opt_output_file, 'r+') as opt_output_file, mmap.mmap(opt_output_file.fileno(), 0, access=mmap.ACCESS_READ) as opt_out:
        if opt_out.find(b'ORCA TERMINATED NORMALLY') != -1:
            print ("ORCA TERMINATED NORMALLY")
            finding=opt_out.find(b'HURRAY')
            if finding==-1:
                print('Geometry NOT optimized')
                with open(opt_input_file,'r+') as opt_file:
                    a=opt_file.readlines()
                    opt_file.seek(0)
                    opt_file.truncate()
                    for line in a:
                        for part in line.split():
                            if ("xyzfile") in part:
                                line=line.strip()
                                line=line.replace(line,"*xyzfile 0 1 "+str(opt_geom_file)+"\n")
                        opt_file.write(line)
                    opt_file.close()    
                with open (opt_error_file,"w") as opt_err:
                    p=sp.Popen(['ORCA',opt_input_file], stdout=opt_output_file, stderr=opt_err)
                    p_status=p.wait()
                    opt_err.close()           
            else:
                print('Geometry optimized successfully')
                copy=False
                for line in opt_output_file:
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
        else:
            print("error")                         
        opt_out.close()
        opt_output_file.close()
        
#generate input file for Freq calculation
freq_keywords_array=[s.replace('Opt' , 'Freq') for s in opt_keywords_array]
freq_geom_array=[s.replace(str(geom_file) , str(opt_geom_file)) for s in geom_array]

with open(freq_input_file, "w") as freq_file:
    for item in freq_keywords_array:
        freq_file.writelines(["%s " %item])
    for item in print_array:
        freq_file.writelines(["%s \n" %item])
    for item in freq_geom_array:
        freq_file.writelines(["%s " %item])
freq_file.close()  

#run Freq calc
freq_out=open(freq_output_file, "w") 
freq_err=open(freq_error_file,"w") 
p=sp.Popen(['ORCA',freq_input_file], stdout=freq_out, stderr=freq_err)
p_status=p.wait()
freq_out.close()
freq_err.close()

#check Freq calc
with open(freq_output_file, 'rb', 0) as freq_out, \
mmap.mmap(freq_out.fileno(), 0, access=mmap.ACCESS_READ) as s:
    if s.find(b'imaginary mode') != -1:
        print('Optimized geom not at minimum')
        with open(opt_geom_file) as f:
            lines=f.readlines()
            lines=lines[2:]
            with open(opt_2_file, "w") as f1:
                f1.writelines(lines)
            f1.close()
            f.close()
    else:
        print('Optimized geom is at global minimum')
freq_out.close()

#generate TDDFT input file
for i in orbital_energies_array [1:]:
    orbital_number=np.array((orbital_energies_array[1:,0]).astype(np.float))
    orbital_energy=np.array((orbital_energies_array[1:,3]).astype(np.float))

for j in edge_data_array:
    element=np.array(edge_data_array[:,0])
    edge=np.array(edge_data_array[:,2])
    energy_theoretical=np.array((edge_data_array[:,4]).astype(np.float))

energy_theoretical_min=energy_theoretical-30
energy_theoretical_max=energy_theoretical

for k in range(0,len(edge)):
    if edge[k]=='K':
        for l in range(0,len(orbital_energy)):
            if -orbital_energy[l] >= energy_theoretical_min[k] and -orbital_energy[l] <= energy_theoretical_max[k]:
                sub_array=[]
                sub_array.append(element[k])
                sub_array.append(orbital_number[l])
                orbital_window_array.append(sub_array)
                
index=np.argwhere(opt_keywords_array=='Opt')
tddft_keywords_array=np.delete(opt_keywords_array,index)
tddft_calc_array=np.array(["\n%tddft","orbwin[0]=0,1,-1,-1","nroots=20","maxdim=200","end"])

with open(tddft_input_file, "w") as tddft_input:
    for item in tddft_keywords_array:
        tddft_input.writelines(["%s " %item])
    for item in print_array:
        tddft_input.writelines(["%s \n" %item])
    for item in tddft_calc_array:
        tddft_input.writelines(["%s \n" %item])
    for item in freq_geom_array:
        tddft_input.writelines(["%s " %item])
tddft_input.close()

#run tddft calc
tddft_out=open(tddft_output_file, "w") 
tddft_err=open(tddft_error_file,"w") 
p=sp.Popen([ORCA,tddft_input_file], stdout=tddft_out, stderr=tddft_err)
p_status=p.wait()
tddft_out.close()
tddft_err.close()
sp.Popen(['orca_mapspc',tddft_output_file,'ABS','-eV','-x0380','-x1410','-n500','-w0.6'])
p_status=p.wait()

'''

stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 
