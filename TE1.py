# -*- coding: utf-8 -*-
"""
Updated on on Fri May 18 11:12:11 2018
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
from collections import defaultdict
import argparse
import sys
#import scriptexit



parser = argparse.ArgumentParser(description='TE1: Theoretical electron density function calulation')

parser.add_argument('in_geom_file',
    type=argparse.FileType('r'),
    help="molecular geometry file to be read in",
    default=sys.stdin, metavar="FILE")

parser.add_argument('-op',
    default="gas",
    choices=['gas', 'solution'],
    help="Select the default set of orca parameters for particular chemical states.")

parser.add_argument("-opi", 
                    dest="file_orca_params", 
                    required=False,
                    help="input file with orca optimation parameters. This over write any default orca optimisation parameters.", 
                    metavar="FILE")

parser.add_argument("-orca", 
                    dest="orca_executable", 
                    required=False,
                    help="path to the orca executable; C:\Orca\orca.exe is the default path", 
                    metavar="FILE")


args = parser.parse_args()

file_geom_read=args.in_geom_file.name


#print("~ Molecular geometry file details: {}".format(args.in_geom_file))
print("~ Molecular geometry file: {}".format(file_geom_read))

print("~ Orca parameter set : {}".format(args.op))

print("~ Orca parameter file: {}".format(args.file_orca_params))

print("~ Orca executable: {}".format(args.orca_executable))

print(args.orca_executable)
print(type(args.orca_executable))

args = parser.parse_args()
print("args:",args)


print("\n\n")

print("START")

ORCA=r"C:\Orca\orca.exe"

if args.orca_executable is None:
    print("The default path for orca, C:\Orca\orca.exe, is used.")

if args.orca_executable is not None:
    ORCA=args.orca_executable
    print("This does not use the default path for orca, instead it used this path: ", ORCA)


path, file_geom = os.path.split(file_geom_read)

print(file_geom)

print(path)


working_dir=os.getcwd()

#file_geom=""


#for file_name in os.listdir(working_dir):
if not file_geom.endswith('.xyz'): 
    sys.exit("ERROR; The molecular geometry file does not have the expected .xyz file extension.")
    
               
index_of_dot = file_geom.index(".") 
file_geom_without_extension = file_geom[:index_of_dot]


resultsdir = r"TE1_"+file_geom_without_extension+r"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

print(resultsdir)
#os.makedirs(mydir)

#path_in=working_dir
path_in=path
#path_out=working_dir+r"\\"+resultsdir
path_out=path+r"\\"+resultsdir
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
opt_error_file=path_out+r"\opt_error.txt"
freq_error_file=path_out+r"\freq_error.txt"


orbital_energies_array=np.array([])
orbital_window_array=[]
edge_data_array=np.array([['C','6','K','0','284.2'],['N','7','K','0','409.9']])

#generate input file for Opt calculation
opt_keywords_gas_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
opt_keywords_solution_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","CPCM(Water)","Grid3", "FinalGrid5", "Opt"])

# default case
opt_keywords_array=opt_keywords_gas_array



if "gas" in args.op:
    print("The default gas orca optimisation parameters are used.")

if "solution" in args.op:
    opt_keywords_array=opt_keywords_solution_array
    print("The default solution orca optimisation parameters are used.")
    
#if "None" in args.file_orca_params:
if args.file_orca_params is None:    
    print("The default orca optimisation parameters are used from the -op flag.")

opt_keywords_infile_array = ["!"]
if args.file_orca_params is not None:    
    array=[]
    filename_desc=args.file_orca_params
    with open(filename_desc) as f_in:
        for line in f_in:
            opt_keywords_infile_array.append(line.strip())
            array.append(line.strip())
    f_in.close()
    opt_keywords_array=opt_keywords_infile_array
    print("The default orca optimisation parameters from the -op flag are not used. The parameters from the orca parameters file are used and these are: ",array)
#print(opt_keywords_infile_array)

#opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
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
    with open(opt_output_file, 'r+') as opt_out_file, mmap.mmap(opt_out_file.fileno(), 0, access=mmap.ACCESS_READ) as opt_out:
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
                opt_err=open(opt_error_file,"w") 
                p=sp.Popen(['ORCA',opt_input_file], stdout=opt_out_file, stderr=opt_err)
                p_status=p.wait()
                opt_err.close()           
            else:
                print('Geometry optimized successfully')
                copy=False
                for line in opt_out_file:
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
        opt_out_file.close()
        
        
        
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
finding=-1
while (finding==-1):
    with open(freq_output_file, 'r') as freq_out_file, mmap.mmap(freq_out_file.fileno(), 0, access=mmap.ACCESS_READ) as freq_out:
        if freq_out.find(b'ORCA TERMINATED NORMALLY') != -1:
            finding=freq_out.find(b'imaginary mode')
            if finding != -1:
                print ('Optimized geom not at minimum')
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
                        #opt_file.close()    
                with open (opt_error_file,"w") as opt_err, open (opt_output_file, "w") as opt_out:
                    p=sp.Popen(['ORCA',opt_input_file], stdout=opt_out, stderr=opt_err)
                    p_status=p.wait()
                    opt_err.close()
                    opt_out.close()
                with open(freq_input_file,'r+') as freq_file:
                    a=freq_file.readlines()
                    freq_file.seek(0)
                    freq_file.truncate()
                    for line in a:
                        for part in line.split():
                            if ("xyzfile") in part:
                                line=line.strip()
                                line=line.replace(line,"*xyzfile 0 1 "+str(opt_geom_file)+"\n")
                        freq_file.write(line)
                with open(freq_error_file,"w") as freq_err:
                    p=sp.Popen(['ORCA',freq_input_file], stdout=freq_out_file, stderr=freq_err)
                    p_status=p.wait()
                    freq_err.close()
            elif finding == -1:
                print('Optimized geom is at global minimum')
                break
        else:
            print ("error")                        
        freq_out.close()
        freq_out_file.close()


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
                sub_array.append(int(orbital_number[l]))
                orbital_window_array.append(sub_array)
d=defaultdict(list)
for lis in orbital_window_array:
    d[lis[0]].append(lis[1],)
orbital_window_array=[list(x for y in i for x in y) for i in d.items()]
                      
index=np.argwhere(opt_keywords_array=='Opt')
tddft_keywords_array=np.delete(opt_keywords_array,index)
for m in range(0,len(orbital_window_array)):
    tddft_input_file=path_out+r"\TDDFT_%s-edge.inp" %orbital_window_array[m][0]
    tddft_output_file=path_out+r"\TDDFT_%s-edge.out" %orbital_window_array[m][0]    
    tddft_error_file=path_out+r"\TDDFT_%s-edge_error.txt" %orbital_window_array[m][0]
    tddft_orbwin_string="orbwin[0]="        
    orbitals=[]
    orbitals.append(orbital_window_array[m][1])
    orbitals.append(orbital_window_array[m][-1])
    tddft_orbwin_string+="%s,%s" %(str(orbitals[0]),str(orbitals[1])) 
    tddft_orbwin_string+=",-1,-1"
    tddft_calc_array=["\n%tddft",tddft_orbwin_string,"nroots=20","maxdim=200","end"]
    with open(tddft_input_file, "w") as tddft_input:
        tddft_input.writelines("#This is %s K-edge calculation \n" %orbital_window_array[m][0])
        for item in tddft_keywords_array:
            tddft_input.writelines(["%s " %item])
        for item in print_array:
            tddft_input.writelines(["%s \n" %item])
        for item in tddft_calc_array:
            tddft_input.writelines(["%s \n" %item])
        for item in freq_geom_array:
            tddft_input.writelines(["%s " %item])
        tddft_input.close()
    #run tddft calculations
    tddft_out=open(tddft_output_file, "w") 
    tddft_err=open(tddft_error_file,"w") 
    p=sp.Popen([ORCA,tddft_input_file], stdout=tddft_out, stderr=tddft_err)
    p_status=p.wait()
    tddft_out.close()
    tddft_err.close()
   
                     
stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 