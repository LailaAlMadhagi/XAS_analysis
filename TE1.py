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


args = parser.parse_args()

file_geom_read=args.in_geom_file.name


#print("~ Molecular geometry file details: {}".format(args.in_geom_file))
print("~ Molecular geometry file: {}".format(file_geom_read))

print("~ Orca parameter set : {}".format(args.op))

print("~ Orca parameter file: {}".format(args.file_orca_params))


print("\n\n")

print("START")

ORCA=r"C:\Orca\orca.exe"


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
tddft_input_file=path_out+r"\TDDFT.inp"
tddft_output_file=path_out+r"\TDDFT.out"
opt_error_file=path_out+r"\opt_error.txt"
freq_error_file=path_out+r"\freq_error.txt"
tddft_error_file=path_out+r"\TDDFT_error.txt"

orbital_energies_array=np.array([])


#generate input file for Opt calculation
opt_keywords_gas_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
opt_keywords_solution_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])

# default case
opt_keywords_array=opt_keywords_gas_array



if "gas" in args.op:
    print("The default gas orca optimisation parameters are used.")

if "solution" in args.op:
    opt_keywords_array=opt_keywords_solution_array
    print("The default solution orca optimisation parameters are used.")
    
#if "None" in args.file_orca_params:
#    print("default orca optimisation parameters are used")

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



stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 
