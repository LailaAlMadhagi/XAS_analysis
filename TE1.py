# -*- coding: utf-8 -*-
"""
Created on Tue May 15 17:51:56 2018
Written by: Laila Al-Madhagi
fy11lham@leeds.ac.uk
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

#resultsdir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
resultsdir = r"results_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

print(resultsdir)
#os.makedirs(mydir)


#working_dir=r"C:\Users\menjle\Desktop\python3_anaconda\PortableGit\NEXAFS_code\gas"
working_dir=os.getcwd()

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
orbital_energies_file=path_out+r"\orbital_energies.txt"
freq_error_file=path_out+r"\freq_error.txt"
tddft_error_file=path_out+r"\TDDFT_error.txt"


#generate input file for Opt calculation
opt_file=open(opt_input_file,"w")
opt_file.write("!")
opt_file.write(" B3LYP")
opt_file.write(" 6-31G*")
opt_file.write(" TIGHTSCF")
opt_file.write(" Grid3 FinalGrid5")
opt_file.write(" Opt ")
opt_file.write("\n")
#print_block="!NormalPrint" "\ %output \n Print[P_Basis] 2 \n Print[P_MOS] 1 \n"
line1="!NormalPrint"
line2="%output"
line3="Print[P_Basis] 2"
line4="Print[P_MOS] 1"
line5="end"
print_block="%s \n%s \n%s \n%s \n%s \n" % (line1, line2, line3, line4,line5)
opt_file.write(print_block)
opt_file.write("*xyzfile 0 1 "+str(geom_file)+"\n")
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
                with open (orbital_energies_file,"w") as orb_eng:
                    copy=False
                    for line in opt_output_file:
                        if line.strip()=="ORBITAL ENERGIES":
                            copy=True
                        elif line.strip()=="MOLECULAR ORBITALS":
                            copy=False
                        elif copy:
                            orb_eng.write(line)
                orb_eng.close()
                with open (orbital_energies_file,"r+") as orb_eng:
                    columns=[]
                    for row in orb_eng:
                        a=row.split()
                        if len(a)==4:
                            columns.append(a)
                    occ_columns=[]
                    x=int(len(columns)/2)
                    for x in range(int(len(columns)/2),len(columns)-1): #this is not true, I need to find a better way to choose the final orbiral energy
                        if columns[x][1]=="0.0000":
                            pass
                        else:
                            occ_columns.append(columns[x])
                        x+=1
                    orb_eng.seek(0)
                    orb_eng.truncate()
                    orb_eng.write("%s \n" %occ_columns) 
                orb_eng.close() 
                
        else:
            print("error")                         
        opt_out.close()
        opt_output_file.close()
#generate input file for Freq calculation
with open(opt_input_file) as opt_file, open(freq_input_file, 'w') as Freq_file:
    for line in opt_file:
        line=line.replace("Opt","Freq")
        for part in line.split():
            if ("xyzfile") in part:
                line=line.strip()
                line=line.replace(line,"*xyzfile 0 1 "+str(opt_geom_file)+"\n")
        Freq_file.write(line)
    opt_file.close()
    Freq_file.close()

#run Freq calc
Freq_out=open(freq_output_file, "w") 
Freq_err=open(freq_error_file,"w") 
p=sp.Popen(['ORCA',freq_input_file], stdout=Freq_out, stderr=Freq_err)
p_status=p.wait()
Freq_out.close()
Freq_err.close()
#check Freq calc
with open(freq_output_file, 'rb', 0) as freq_output_file, \
     mmap.mmap(freq_output_file.fileno(), 0, access=mmap.ACCESS_READ) as s:
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
freq_output_file.close()

#generate TDDFT input file
with open(freq_input_file,"r") as Freq_file, open(tddft_input_file, 'w') as tddft_file:
    for line in Freq_file:
        line=line.replace("Freq","")
        tddft_file.write(line)
with open(tddft_input_file,"r+") as tddft_file:
    a=tddft_file.readlines()
    index=0
    for item in a:
        if item.startswith("!Normal"):
            tddft_line1="\n%tddft"
            tddft_line2="orbwin[0]=0,1,-1,-1"
            tddft_line3="nroots=20"
            tddft_line4="maxdim=200"
            tddft_line5="end"
            a.insert(index,"%s \n%s \n%s \n%s \n%s \n" % (tddft_line1, tddft_line2, tddft_line3, tddft_line4, tddft_line5))
            break
        index+=1
    tddft_file.seek(0)
    tddft_file.truncate()
    for line in a:
        tddft_file.write(line)
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