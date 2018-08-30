# -*- coding: utf-8 -*-
"""
Updated on on Fri August 20 11:17:11 2018
Written by: Laila Al-Madhagi on May 11 2018
fy11lham@leeds.ac.uk
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
description='LES1: Liquid Excited State calulation'
parser = argparse.ArgumentParser(description)


parser.add_argument('in_geom_dir',
    #type=argparse.FileType('r'),
    help="molecular geometry file to be read in",
    default=sys.stdin, metavar="DIR")

parser.add_argument("-geom_file_name",
                    dest="geom_file_name",
                    type=str,
                    required=False,
                    help="the geometry file name")
"""
parser.add_argument('-op',
    default="solution",
    help="Select the default set of orca parameters for particular chemical states.")
"""

parser.add_argument("-opi", 
                    dest="file_orca_params", 
                    required=False,
                    help="input file with orca parameters. This over write any default orca parameters.", 
                    metavar="FILE")

parser.add_argument("-orca", 
                    dest="orca_executable", 
                    required=False,
                    help="path to the orca executable; C:\Orca\orca.exe is the default path", 
                    metavar="FILE")

parser.add_argument("-path_out",
                    dest="LES_path_out",
                    type=str,
                    required=False,
                    help="directory where output is written")

parser.add_argument("-element",
                    dest="element",
                    type=str,
                    required=False,
                    help="the elemtent for which excited state calculations will be performed")

parser.add_argument("-pal",
                    dest="processors_number",
                    type=str,
                    required=False,
                    help="If running ORCA in parallel specify the number of processors")


# all input argument have been read now we can process them
args = parser.parse_args()

if args.in_geom_dir is None:
    if os.path.split(args.geom_file_name)[0] == '':
        sys.exit("ERROR; There is no geometry file provided")
    elif os.path.split(args.geom_file_name)[0] != '':
        geom_files_dir=os.path.split(args.geom_file_name)[0]

path_LES, geom_files_dir = os.path.split(args.in_geom_dir)
LES_date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

working_dir=os.getcwd()

if args.geom_file_name is not None:
    resultsdir = r"LES1_"+args.geom_file_name+r"_"+LES_date_time
    geom_file=args.in_geom_dir+r'//'+args.geom_file_name
if args.geom_file_name is None:
    geom_file=args.in_geom_dir+r'//'+os.listdir(args.in_geom_dir)[0]
    resultsdir = r"LES1_"+geom_file.split('//')[1]+r'_'+LES_date_time
  
# set in and out paths   
script_path=os.path.dirname(os.path.abspath(__file__))    
path_in=working_dir
if args.LES_path_out is not None:
    path_out=args.LES_path_out+r"//"+resultsdir
if args.LES_path_out is None:
    path_out=path_in+r"//"+resultsdir
os.makedirs(path_out)

# set log file
log_file_name = path_out+r"//log.txt"
log_file=open(log_file_name, "w") 
log_file.write(description+"\n\n")
try:
    host=socket.gethostbyaddr(socket.gethostname())[0]
except socket.herror:
    host=''
log_file.write(r"This program ran at "+SES_date_time+r" on the "+host+r" host system.")

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


print("\n~ Molecular geometry file details: {}".format(args.in_geom_dir))
log_file.write("\n\n~ Molecular geometry file details: {}".format(args.in_geom_dir))

"""
print("\n~ Orca parameter set : {}".format(args.op))
log_file.write("\n\n~ Orca parameter set : {}".format(args.op))
"""

print("\n~ Orca parameter file: {}".format(args.file_orca_params))
log_file.write("\n\n~ Orca parameter file: {}".format(args.file_orca_params))

ORCA=r"C://Orca//orca"

if args.orca_executable is None or args.orca_executable=='':
    print("The default path for orca, C://Orca//orca, is used.")
    log_file.write("\n\nThe default path for orca, C://Orca//orca, is used.")

if args.orca_executable is not None and args.orca_executable!='':
    ORCA=args.orca_executable
    print("This does not use the default path for orca, instead it used this path: "+ORCA)
    log_file.write("\n\nThis does not use the default path for orca, instead it used this path: "+ORCA)



print("\n\nSTART: \n")
log_file.write("\n\nSTART: \n")
log_file.flush()


if args.processors_number is None or args.processors_number == '':
    parallel_array=np.array(["#This is a serial job","\n\n"])      
if args.processors_number is not None and args.processors_number != '':
    parallel_array=np.array(["#This is a parallel job", "\n%pal"," nprocs ",int(args.processors_number)," end"])

    
# Here are all the files we need to create (except the log file which we are already using)

sp_input_file=path_out+r"//sp.inp"
new_sp_input_file=path_out+r"//new_sp.inp"
sp_output_file=path_out+r"//sp.out"
sp_error_file=path_out+r"//sp_error.txt"
edge_data_file=script_path+r"//..//edge_data.txt"
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

#generate input file for Single point calculation
# default case for ORCA SP parameters
sp_keywords_solution_array=np.array(["\n!","B3LYP","D3","ma-def2-SVP","TIGHTSCF","Grid3", "FinalGrid5"])

"""
if "solution" in args.op:
    sp_keywords_array=sp_keywords_solution_array
    print("The default solution orca parameters are used.")
    log_file.write("\n\nThe default solution orca parameters are used.")
"""

if args.file_orca_params is None:
    sp_keywords_array=sp_keywords_solution_array    
    print("The default orca optimisation parameters are used from the -op flag.")
    log_file.write("\n\nThe default orca optimisation parameters are used from the -op flag.")
log_file.flush()

#user input for ORCA SP parameters
sp_keywords_infile_array = ["\n!"]
if args.file_orca_params is not None:    
    array=[]
    filename_desc=args.file_orca_params
    with open(filename_desc) as f_in:
        for line in f_in:
            sp_keywords_infile_array.append(line.strip())
            array.append(line.strip())
    f_in.close()
    sp_keywords_array=sp_keywords_infile_array
    print("The default orca optimisation parameters from the -op flag are not used. Instead the parameters from the orca parameters file are used and these are: ")
    print(array)
    log_file.write("\n\nThe default orca optimisation parameters from the -op flag are not used. Instead the parameters from the orca parameters file are used and these are: \n")
    log_file.write("\n".join(str(elem) for elem in array))

#ORCA parameters for printing block
print_array=np.array(["\n!NormalPrint","%output","Print[P_Basis] 2","Print[P_MOS] 1","end"])
print("\n")
log_file.write("\n\n")

#ORCA parameters for geometry block
geom_array=np.array(["*xyzfile","0","1",str(geom_file)])


#write SP input file
with open(sp_input_file, "w") as sp_file:
    for item in parallel_array:
        sp_file.writelines(["%s " %item])
    for item in sp_keywords_array:
        sp_file.writelines(["%s " %item])
    for item in print_array:
        sp_file.writelines(["%s \n" %item])
    for item in geom_array:
        sp_file.writelines(["%s " %item])
sp_file.close()

# run sp calculation
sp_out=open(sp_output_file, "w") 
sp_err=open(sp_error_file,"w") 
try:
    sp_p=sp.Popen([str(ORCA),sp_input_file], stdout=sp_out, stderr=sp_err)
    sp_p_status=sp_p.wait()
    sp_out.close()
    sp_err.close()
except Exception as e:
    exc_type,exc_obj,exc_tb=sys.exc_info()
    line = exc_tb.tb_lineno
    fname=exc_tb.tb_next.tb_frame.f_code.co_filename
    sys.exit("Error calling ORCA using subprocess. Error type is {0}, check line {1} in code".format(exc_type,line))


log_file.close()
log_file=open(log_file_name, "w")

# check sp.out file
with open(sp_output_file, 'r+') as sp_out_file, mmap.mmap(sp_out_file.fileno(), 0, access=mmap.ACCESS_READ) as sp_out:
    if sp_out.find(b'ORCA TERMINATED NORMALLY') != -1:
        print ("ORCA TERMINATED NORMALLY")
        log_file.write("ORCA TERMINATED NORMALLY\n")
        copy=False
        for line in sp_out_file:
            if line.strip()=="ORBITAL ENERGIES":
                copy=True
            elif line.strip()=="MOLECULAR ORBITALS":
                copy=False
            elif copy:
                line=line.split()
                if len(line)==4:
                    orbital_energies_array=np.append(orbital_energies_array,line)
        orbital_energies_array=orbital_energies_array.reshape(int(len(orbital_energies_array)/4),4)
        unocc_orbitals=[]
        for x in range (0, int(len(orbital_energies_array)-1)):
            if orbital_energies_array[x][1]==np.str("0.0000"):
                unocc_orbitals.append(x)
            x+=1
        orbital_energies_array=np.delete(orbital_energies_array, np.s_[unocc_orbitals[0]:],0)            
    else:
        print("Error in single point energy calculations.") 
        log_file.write("single point energy calculations.\n")  
    sp_out.close()
    sp_out_file.close()
print("\nEnd of Single Point Energy Calculation.\n\n ********************** \n")
log_file.write("\nEnd of Single Point Energy Calculation.\n\n ********************** \n")
log_file.flush()

sp_timer = timeit.default_timer() 

    
#list of elements from geom file
elements_geom_in=[]
with open (geom_file,'r') as geom_in:
    for line in geom_in:
        elements_geom_in.append(line[0])


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
tddft_keywords_array=sp_keywords_array
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
        for item in parallel_array:
            tddft_input.writelines(["%s " %item])
        for item in tddft_keywords_array:
            tddft_input.writelines(["%s " %item])
        for item in print_array:
            tddft_input.writelines(["%s \n" %item])
        for item in tddft_calc_array:
            tddft_input.writelines(["%s \n" %item])
        for item in geom_array:
            tddft_input.writelines(["%s " %item])
        tddft_input.close()
    #run tddft calculations
    tddft_out=open(tddft_output_file, "w") 
    tddft_err=open(tddft_error_file,"w") 
    try:
        tddft_p=sp.Popen([str(ORCA),tddft_input_file], stdout=tddft_out, stderr=tddft_err)
        tddft_p_status=tddft_p.wait()
        tddft_out.close()
        tddft_err.close()
    except Exception as e:
        exc_type,exc_obj,exc_tb=sys.exc_info()
        line = exc_tb.tb_lineno
        fname=exc_tb.tb_next.tb_frame.f_code.co_filename
        sys.exit("Error calling ORCA using subprocess. Error type is {0}, check line {1} in code".format(exc_type,line))

log_file.close()
log_file=open(log_file_name, "w")
log_file.flush()
   
#print running time information                   
stop = timeit.default_timer()
running_time=(stop-start)/60
print ("\n\nEND: \nRunning time is: "+ str(round(running_time,3)) + " minutes") 
log_file.write("\n\nEND:\nRunning time is: "+ str(round(running_time,3)) + " minutes")
sp_time=(sp_timer-start)/60
tddft_time=(stop-sp_timer)/60
print("\tSingle Point time is: "+ str(round(sp_time,3)) + " minutes")
print("\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")
log_file.write("\n\tSingle Point Energy time is: "+ str(round(sp_time,3)) + " minutes")
log_file.write("\n\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")

print ("LES1 Process Ended Successfully")
log_file.write("LES1 Process Ended Successfully")
print("\n~ path for directory where outputs are: {}".format(path_out))
log_file.write("\n~ path for directory where outputs are: {}".format(path_out))


log_file.close()
            