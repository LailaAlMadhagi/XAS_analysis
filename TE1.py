# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:42:12 2018

Written by: Laila Al-Madhagi
fy11lham@leeds.ac.uk
"""
import timeit
start = timeit.default_timer()
import subprocess as sp
import mmap
import numpy as np

ORCA=r"C:\Orca\orca.exe"
working_dir=r"C:\Users\dmg81179\Desktop\Code_Development\2018-March-28_ORCA_code\working_dir"
geom_file=working_dir+r"\geom.xyz"
Opt_input_file=working_dir+r"\Opt.inp"
new_Opt_input_file=working_dir+r"\new_Opt.inp"
Opt_output_file=working_dir+r"\Opt.out"
new_Opt_output_file=working_dir+r"\new_Opt.out"
Opt_geom_file=working_dir+r"\Opt.xyz"
Freq_input_file=working_dir+r"\Freq.inp"
Freq_output_file=working_dir+r"\Freq.out"
tddft_input_file=working_dir+r"\tddft.inp"
tddft_output_file=working_dir+r"\tddft.out"
orbital_energies_array=np.array([])

#generate input file for Opt calculation
opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
print_array=np.array(["\n","!NormalPrint","%output","Print[P_Basis] 2","Print[P_MOS] 1","end"])
geom_array=np.array(["*xyzfile","0","1",str(geom_file)])

with open(Opt_input_file, "w") as Opt_file:
    for item in opt_keywords_array:
        Opt_file.writelines(["%s " %item])
    for item in print_array:
        Opt_file.writelines(["%s \n" %item])
    for item in geom_array:
        Opt_file.writelines(["%s " %item])
Opt_file.close()

# run Opt calculation
Opt_out=open(Opt_output_file, "w") 
Opt_err=open(working_dir+"\Opt_error.txt","w") 
p=sp.Popen(['ORCA',Opt_input_file], stdout=Opt_out, stderr=Opt_err)
p_status=p.wait()
Opt_out.close()
Opt_err.close()

# check Opt.out file
finding=-1
while (finding==-1):
    with open(Opt_output_file, 'r+') as Opt_output_file, mmap.mmap(Opt_output_file.fileno(), 0, access=mmap.ACCESS_READ) as Opt_out:
        if Opt_out.find(b'ORCA TERMINATED NORMALLY') != -1:
            print ("ORCA TERMINATED NORMALLY")
            finding=Opt_out.find(b'HURRAY')
            if finding==-1:
                print('Geometry NOT optimized')
                with open(Opt_input_file,'r+') as Opt_file:
                    a=Opt_file.readlines()
                    Opt_file.seek(0)
                    Opt_file.truncate()
                    for line in a:
                        for part in line.split():
                            if ("xyzfile") in part:
                                line=line.strip()
                                line=line.replace(line,"*xyzfile 0 1 "+str(Opt_geom_file)+"\n")
                        Opt_file.write(line)
                    Opt_file.close()    
                with open (working_dir+r"\Opt_error.txt","w") as Opt_err:
                    p=sp.Popen(['ORCA',Opt_input_file], stdout=Opt_output_file, stderr=Opt_err)
                    p_status=p.wait()
                    Opt_err.close()
            
            else:
                copy=False
                for line in Opt_output_file:
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
                print('Geometry optimized successfully')
                
        else:
            print("error")                         
        Opt_out.close()
        Opt_output_file.close()
        
#generate input file for Freq calculation
freq_keywords_array=[s.replace('Opt' , 'Freq') for s in opt_keywords_array]
freq_geom_array=[s.replace(str(geom_file) , str(Opt_geom_file)) for s in geom_array]

with open(Freq_input_file, "w") as freq_file:
    for item in freq_keywords_array:
        freq_file.writelines(["%s " %item])
    for item in print_array:
        freq_file.writelines(["%s \n" %item])
    for item in freq_geom_array:
        freq_file.writelines(["%s " %item])
freq_file.close()  

#run Freq calc
Freq_out=open(Freq_output_file, "w") 
Freq_err=open(working_dir+"\Freq_error.txt","w") 
p=sp.Popen(['ORCA',Freq_input_file], stdout=Freq_out, stderr=Freq_err)
p_status=p.wait()
Freq_out.close()
Freq_err.close()
#check Freq calc
with open(Freq_output_file, 'rb', 0) as Freq_output_file, \
     mmap.mmap(Freq_output_file.fileno(), 0, access=mmap.ACCESS_READ) as s:
    if s.find(b'imaginary mode') != -1:
        print('Optimized geom not at minimum')
        with open(Opt_geom_file) as f:
            lines=f.readlines()
            lines=lines[2:]
            with open(working_dir+"\Opt_2.xyz", "w") as f1:
                f1.writelines(lines)
            f1.close()
            f.close()
    else:
        print('Optimized geom is at global minimum')
Freq_output_file.close()

#generate TDDFT input file
index=np.argwhere(opt_keywords_array=='Opt')
TDDFT_keywords_array=np.delete(opt_keywords_array,index)
TDDFT_calc_array=np.array(["\n%tddft","orbwin[0]=0,1,-1,-1","nroots=20","maxdim=200","end"])

with open(tddft_input_file, "w") as TDDFT_file:
    for item in TDDFT_keywords_array:
        TDDFT_file.writelines(["%s " %item])
    for item in print_array:
        TDDFT_file.writelines(["%s \n" %item])
    for item in TDDFT_calc_array:
        TDDFT_file.writelines(["%s \n" %item])
    for item in freq_geom_array:
        TDDFT_file.writelines(["%s " %item])
TDDFT_file.close()

#run tddft calc
tddft_out=open(tddft_output_file, "w") 
tddft_err=open(working_dir+"\error_tddft.txt","w") 
p=sp.Popen([ORCA,tddft_input_file], stdout=tddft_out, stderr=tddft_err)
p_status=p.wait()
tddft_out.close()
tddft_err.close()
sp.Popen(['orca_mapspc',tddft_output_file,'ABS','-eV','-x0380','-x1410','-n500','-w0.6'])
p_status=p.wait()


stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 