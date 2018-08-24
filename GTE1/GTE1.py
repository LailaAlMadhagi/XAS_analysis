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
import socket




# handle the input flags
description='GTE1: Theoretical electron density function calulation for gas'
parser = argparse.ArgumentParser(description)


parser.add_argument('in_geom_file',
    type=argparse.FileType('r'),
    help="molecular geometry file to be read in",
    default=sys.stdin, metavar="FILE")

'''parser.add_argument("-v", 
                    "--verbose", 
                    const=1, 
                    default=0, 
                    type=int, 
                    nargs="?",
                    help="increase verbosity: 0 = only warnings, 1 = info, 2 = debug, 3 = write to log file in results directory. No number means info. Default is no verbosity.")
'''

parser.add_argument("-opi", 
                    dest="file_orca_params", 
                    required=False,
                    help="input file with orca optimation parameters. This over write any default orca optimisation parameters.", 
                    metavar="FILE")

parser.add_argument("-orca", 
                    dest="orca_executable",
                    type=str,
                    required=False,
                    help="path to the orca executable; C:\Orca\orca.exe is the default path")

parser.add_argument("-path_out",
                    dest="path_out",
                    type=str,
                    required=False,
                    help="directory where output is written")

parser.add_argument("-element",
                    dest="element",
                    type=str,
                    required=False,
                    help="the elemtent for which excited state calculations will be performed")

# all input argument have been read now we can process them
args = parser.parse_args()

file_geom_read=args.in_geom_file.name

path, file_geom = os.path.split(file_geom_read)



if not file_geom.endswith('.xyz'): 
    sys.exit("ERROR; The molecular geometry file does not have the expected .xyz file extension.")
                   
index_of_dot = file_geom.index(".") 
file_geom_without_extension = file_geom[:index_of_dot]

date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
resultsdir = r"GTE1_"+file_geom_without_extension+r"_"+date_time


path_in=path
if args.path_out is not None:
    path_out=args.path_out+r"//"+resultsdir
if args.path_out is None:
    path_out=path_in+r"//"+resultsdir
os.makedirs(path_out)

# We do not use the logfile handler because it conflicts with the argparser that
# controls the command line options.
log_file_name = path_out+r"//%s_GTE1_log.txt"%file_geom_without_extension
log_file=open(log_file_name, "w") 
log_file.write(description+"\n\n")

try:
    host=socket.gethostbyaddr(socket.gethostname())[0]
except socket.herror:
    host=''
log_file.write(r"This program ran at "+date_time+r" on the "+host+r" host system.")

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

log_file.write("\n System's IP address is: "+IP)
log_file.write("\n\n")


print("\n~ Molecular geometry file details: {}".format(args.in_geom_file))
log_file.write("\n\n~ Molecular geometry file details: {}".format(args.in_geom_file))

print("\n~ Orca parameter file: {}".format(args.file_orca_params))
log_file.write("\n\n~ Orca parameter file: {}".format(args.file_orca_params))

print("\n~ Orca executable: {}".format(args.orca_executable))
log_file.write("\n\n~ Orca executable: {}".format(args.orca_executable))

print("\n~ Excited state calculations are run for element: {}".format(args.element))
log_file.write("\n~ Excited state calculations are run for element: {}".format(args.element))

log_file.write('\n\nPython {0} and {1}'.format((sys.version).split('|')[0],(sys.version).split('|')[1]))
log_file.write('\n\nNumpy version is: '+np.__version__)
log_file.flush()

#print(args.orca_executable)
#print(type(args.orca_executable))

#args = parser.parse_args()
#args_text = "args:".join(args)


print("\n\nSTART: \n")
log_file.write("\n\nSTART: \n")
log_file.flush()

ORCA=r"C://Orca//orca"

if args.orca_executable is None or args.orca_executable=='':
    print("The default path for orca, C://Orca//orca, is used.")
    log_file.write("\n\nThe default path for orca, C://Orca//orca, is used.")

if args.orca_executable is not None and args.orca_executable!='':
    ORCA=args.orca_executable
    print("This does not use the default path for orca, instead it used this path: "+ORCA)
    log_file.write("\n\nThis does not use the default path for orca, instead it used this path: "+ORCA)



# Here are all the files we need to create (except the log file which we are already using)
geom_file=args.in_geom_file.name
opt_input_file=path_out+r"//opt.inp"
new_opt_input_file=path_out+r"//new_opt.inp"
opt_output_file=path_out+r"//opt.out"
new_opt_output_file=path_out+r"//new_Opt.out"
opt_geom_file=path_out+r"//opt.xyz"
opt_2_file=path_out+r"//opt_2.xyz"
freq_input_file=path_out+r"//freq.inp"
freq_output_file=path_out+r"//freq.out"
opt_error_file=path_out+r"//opt_error.txt"
freq_error_file=path_out+r"//freq_error.txt"
working_dir=os.getcwd()
edge_data_file=working_dir+r"//..//edge_data.txt"

orbital_energies_array=np.array([])
all_orbital_window_array=[]
final_orbital_window_array=[]
#edge_data_array=np.array([['C','6','K','0','284.2'],['N','7','K','0','409.9']])

#put edge data into array
edge_data_array=np.array([])
with open(edge_data_file,"r") as edge_data_file:
    next (edge_data_file)
    next (edge_data_file)
    for line in edge_data_file:
        line=line.split()
        edge_data_array=np.append(edge_data_array,line)
    edge_data_array=edge_data_array.reshape(int(len(edge_data_array)/5),5)
    edge_data_file.close()



#generate input file for Opt calculation
opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])

if args.file_orca_params is None:
    print("The default gas orca optimisation parameters are used.")
    log_file.write("\n\nThe default gas orca optimisation parameters are used.")

log_file.flush()

#list of elements from geom file
elements_geom_in=[]
with open (geom_file,'r') as geom_in:
    for line in geom_in:
        elements_geom_in.append(line[0])

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
    print("The default orca optimisation parameters are not used. Instead the parameters from the orca parameters file are used and these are: ")
    print(array)
    log_file.write("\n\nThe default orca optimisation parameters are not used. Instead the parameters from the orca parameters file are used and these are: \n")
    log_file.write("\n".join(str(elem) for elem in array))

#opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
print_array=np.array(["\n!NormalPrint","%output","Print[P_Basis] 2","Print[P_MOS] 1","end"])
print("\n")
log_file.write("\n\n")
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
p=sp.Popen([str(ORCA),opt_input_file], stdout=opt_out, stderr=opt_err)
p_status=p.wait()
opt_out.close()
opt_err.close()

# The log_file write times out after about an hour. We close it and then
# open them again because this loop can take a long time to run
log_file.close()
log_file=open(log_file_name, "w")
# check Opt.out file
loop=1
finding=-1
while (finding==-1):
    print(r"Loop count is "+str(loop)+r".")
    log_file.write("\n\nLoop count is "+str(loop)+".\n\n")
    loop+=1
    with open(opt_output_file, 'r+') as opt_out_file, mmap.mmap(opt_out_file.fileno(), 0, access=mmap.ACCESS_READ) as opt_out:
        if opt_out.find(b'ORCA TERMINATED NORMALLY') != -1:
            print ("ORCA TERMINATED NORMALLY")
            log_file.write("ORCA TERMINATED NORMALLY\n")
            finding=opt_out.find(b'HURRAY')
            if finding==-1:
                print('Geometry NOT optimized')
                log_file.write('Geometry NOT optimized\n')
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
                p=sp.Popen([str(ORCA),opt_input_file], stdout=opt_out_file, stderr=opt_err)
                p_status=p.wait()
                opt_err.close()           
            else:
                print('Geometry optimized successfully')
                log_file.write('Geometry optimized successfully\n')
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
            print("Error in optimisation loop.") 
            log_file.write("Error in optimisation loop.\n")                        
        opt_out.close()
        opt_out_file.close()
    print("\nEnd of Optimisation loop.\n\n ********************** \n")
    log_file.write("\nEnd of Optimisation loop.\n\n ********************** \n")
    log_file.flush()
        
opt_timer = timeit.default_timer()        
        
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
p=sp.Popen([str(ORCA),freq_input_file], stdout=freq_out, stderr=freq_err)
p_status=p.wait()
freq_out.close()
freq_err.close()
log_file.close()
log_file=open(log_file_name, "w")
#check Freq calc
finding=0
while (finding==0):
    with open(freq_output_file, 'r+') as freq_out_file, mmap.mmap(freq_out_file.fileno(), 0, access=mmap.ACCESS_READ) as freq_out:
        if freq_out.find(b'ORCA TERMINATED NORMALLY') != -1:
            finding=freq_out.find(b'imaginary mode')
            if finding == 0:
                print ('Optimized geom not at minimum')
                log_file.write('Optimized geom not at minimum\n')
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
                with open (opt_error_file,"w") as opt_err, open (opt_output_file, "w") as opt_out:
                    p=sp.Popen([str(ORCA),opt_input_file], stdout=opt_out, stderr=opt_err)
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
                freq_err=open(freq_error_file,"w")
                p=sp.Popen([str(ORCA),freq_input_file], stdout=freq_out_file, stderr=freq_err)
                p_status=p.wait()
                freq_err.close()
            elif finding !=0:
                print('Optimized geom is at global minimum')
                log_file.write('Optimized geom is at global minimum\n')
                break
        else:
            print ("Error: in freq loop")
            log_file.write("Error: in freq loop")                   
        freq_out.close()
        freq_out_file.close()

freq_timer = timeit.default_timer() 
log_file.flush()

#generate TDDFT input file
for i in orbital_energies_array [1:]:
    orbital_number=np.array((orbital_energies_array[1:,0]).astype(np.float))
    orbital_energy=np.array((orbital_energies_array[1:,3]).astype(np.float))

for j in edge_data_array:
    element=np.array(edge_data_array[:,0])
    edge=np.array(edge_data_array[:,2])
    energy_theoretical=np.array((edge_data_array[:,4]).astype(np.float))

energy_theoretical_min=energy_theoretical-60
energy_theoretical_max=energy_theoretical

for k in range(0,len(edge)):
    if edge[k]=='K':
        for l in range(0,len(orbital_energy)):
            if -orbital_energy[l] >= energy_theoretical_min[k] and -orbital_energy[l] <= energy_theoretical_max[k]:
                sub_array=[]
                sub_array.append(element[k])
                sub_array.append(int(orbital_number[l]))
                if sub_array[0] in elements_geom_in:
                    all_orbital_window_array.append(sub_array)
                else:
                    pass
d=defaultdict(list)
for lis in all_orbital_window_array:
    d[lis[0]].append(lis[1],)
all_orbital_window_array=[list(x for y in i for x in y) for i in d.items()]
                      
index=np.argwhere(opt_keywords_array=='Opt')
tddft_keywords_array=np.delete(opt_keywords_array,index)
for i in all_orbital_window_array:
    if args.element is None:
        final_orbital_window_array.append(i)        
    if args.element is not None and i[0]==args.element:     
        final_orbital_window_array.append(i)
for j in final_orbital_window_array:
    tddft_input_file=path_out+r"//TDDFT_%s-edge.inp" %j[0]
    tddft_output_file=path_out+r"//TDDFT_%s-edge.out" %j[0]    
    tddft_error_file=path_out+r"//TDDFT_%s-edge_error.txt" %j[0]
    tddft_orbwin_string="orbwin[0]="        
    orbitals=[]
    orbitals.append(j[1])
    orbitals.append(j[-1])
    tddft_orbwin_string+="%s,%s" %(str(orbitals[0]),str(orbitals[1])) 
    tddft_orbwin_string+=",-1,-1"
    tddft_calc_array=["\n%tddft",tddft_orbwin_string,"nroots=20","maxdim=200","end"]
    with open(tddft_input_file, "w") as tddft_input:
        tddft_input.writelines("#This is %s K-edge calculation \n" %j[0])
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
    p=sp.Popen([str(ORCA),tddft_input_file], stdout=tddft_out, stderr=tddft_err)
    p_status=p.wait()
    tddft_out.close()
    tddft_err.close()

log_file.close()
log_file=open(log_file_name, "w")
log_file.flush()
                     
stop = timeit.default_timer()
running_time=(stop-start)/60
print ("\n\nEND: \nRunning time is: "+ str(round(running_time,3)) + " minutes") 

log_file.write("\n\nEND:\nRunning time is: "+ str(round(running_time,3)) + " minutes")

opt_time=(opt_timer-start)/60
freq_time=(freq_timer-opt_timer)/60
tddft_time=(stop-freq_timer)/60

print("\tOptimisation time is: "+ str(round(opt_time,3)) + " minutes")
print("\tFrequency calculation time is: "+ str(round(freq_time,3)) + " minutes")
print("\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")

log_file.write("\n\tOptimisation time is: "+ str(round(opt_time,3)) + " minutes")
log_file.write("\n\tFrequency calculation time is: "+ str(round(freq_time,3)) + " minutes")
log_file.write("\n\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")

print("\n~ path for directory where outputs are: {}".format(path_out))
log_file.write("\n~ path for directory where outputs are: {}".format(path_out))

log_file.close()