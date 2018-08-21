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
description='LES: Liquid Excited State calulation'
parser = argparse.ArgumentParser(description)


parser.add_argument('in_geom_dir',
    #type=argparse.FileType('r'),
    help="molecular geometry file to be read in",
    default=sys.stdin, metavar="DIR")

parser.add_argument('-op',
    default="solution",
    help="Select the default set of orca parameters for particular chemical states.")

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


# all input argument have been read now we can process them
args = parser.parse_args()
path_LES, file_geom_dir = os.path.split(args.in_geom_dir)
#file_geom_read=args.in_geom_file.name
geom_filenames=[]
for filename in os.listdir(path_LES+r"/%s"%file_geom_dir):
    geom_filenames.append(filename)
"""
#run E2
p=sp.Popen(['python','U:\XAS_analysis\E2\E2.py',
            'U:\XAS_analysis\E2\Imdz_Solution_35C_APS.txt','0','3','7','--file_type', 'user_defined',
            '-offset', '38'])
p_status=p.wait()
"""    
#if not file_geom.endswith('.xyz'): 
#    sys.exit("ERROR; The molecular geometry file does not have the expected .xyz file extension.")
                 
#index_of_dot = file_geom.index(".") 
#file_geom_without_extension = file_geom[:index_of_dot]

for file in geom_filenames:
    date_time = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    resultsdir = r"LES_"+file_geom_dir+r"_"+file+r"_"+date_time
    
    
    path_in_LES=path_LES
    path_out_LES=path_LES+r"/"+resultsdir
    os.makedirs(path_out_LES)
    
    log_file_name = path_out_LES+r"\\log.txt"
    log_file=open(log_file_name, "w") 
    
    log_file.write(description+"\n\n")
    host=socket.gethostbyaddr(socket.gethostname())[0]
    log_file.write(r"This program ran at "+date_time+r" on the "+host+r" host system.")
    log_file.write("\n\n")
    
    
    
    
    print("\n~ Molecular geometry file details: {}".format(args.in_geom_dir))
    log_file.write("\n\n~ Molecular geometry file details: {}".format(args.in_geom_dir))
    
    
    print("\n~ Orca parameter set : {}".format(args.op))
    log_file.write("\n\n~ Orca parameter set : {}".format(args.op))
    
    
    print("\n~ Orca parameter file: {}".format(args.file_orca_params))
    log_file.write("\n\n~ Orca parameter file: {}".format(args.file_orca_params))
    
    print("\n~ Orca executable: {}".format(args.orca_executable))
    log_file.write("\n\n~ Orca executable: {}".format(args.orca_executable))
    
    
    #print(args.orca_executable)
    #print(type(args.orca_executable))
    
    #args = parser.parse_args()
    #args_text = "args:".join(args)
    
    
    print("\n\nSTART: \n")
    log_file.write("\n\nSTART: \n")
    log_file.flush()
    
    ORCA=r"C:\Orca\orca.exe"
    
    if args.orca_executable is None:
        print("The default path for orca, C:\Orca\orca.exe, is used.")
        log_file.write("\n\nThe default path for orca, C:\Orca\orca.exe, is used.")
    
    if args.orca_executable is not None:
        ORCA=args.orca_executable
        print("This does not use the default path for orca, instead it used this path: ", ORCA)
        log_file.write("\n\nThis does not use the default path for orca, instead it used this path: ", ORCA)



# Here are all the files we need to create (except the log file which we are already using)
    geom_file=path_LES+r"/%s"%file_geom_dir+"/%s"%file
    sp_input_file=path_out_LES+r"\sp.inp"
    new_sp_input_file=path_out_LES+r"\new_sp.inp"
    sp_output_file=path_out_LES+r"\sp.out"
    new_sp_output_file=path_out_LES+r"\new_sp.out"
    sp_error_file=path_out_LES+r"\sp_error.txt"
    working_dir=os.getcwd()
    edge_data_file=working_dir+r"\..\edge_data.txt"
    orbital_energies_array=np.array([])
    orbital_window_array=[]
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
    sp_keywords_solution_array=np.array(["!","B3LYP","D3","ma-def2-SVP","TIGHTSCF","Grid3", "FinalGrid5"])
    
    # default case
    sp_keywords_array=sp_keywords_solution_array
    
    
    
    if "solution" in args.op:
        sp_keywords_array=sp_keywords_solution_array
        print("The default solution orca parameters are used.")
        log_file.write("\n\nThe default solution orca parameters are used.")
        
    #if "None" in args.file_orca_params:
    if args.file_orca_params is None:    
        print("The default orca optimisation parameters are used from the -op flag.")
        log_file.write("\n\nThe default orca optimisation parameters are used from the -op flag.")
    log_file.flush()
    
    #list of elements from geom file
    elements_geom_in=[]
    with open (geom_file,'r') as geom_in:
        for line in geom_in:
            elements_geom_in.append(line[0])
    
    sp_keywords_infile_array = ["!"]
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
    
    #opt_keywords_array=np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])np.array(["!","B3LYP","6-31G*","TIGHTSCF","Grid3", "FinalGrid5", "Opt"])
    
    print_array=np.array(["\n!NormalPrint","%output","Print[P_Basis] 2","Print[P_MOS] 1","end"])
    print("\n")
    log_file.write("\n\n")
    geom_array=np.array(["*xyzfile","0","1",str(geom_file)])
    
    with open(sp_input_file, "w") as sp_file:
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
    p=sp.Popen(['ORCA',sp_input_file], stdout=sp_out, stderr=sp_err)
    p_status=p.wait()
    sp_out.close()
    sp_err.close()
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
            print("Error in optimisation loop.") 
            log_file.write("Error in optimisation loop.\n")  
        sp_out.close()
        sp_out_file.close()
    print("\nEnd of Single Point Energy Calculation.\n\n ********************** \n")
    log_file.write("\nEnd of Single Point Energy Calculation.\n\n ********************** \n")
    log_file.flush()
            
    sp_timer = timeit.default_timer()        
    
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
                        orbital_window_array.append(sub_array)
                    else:
                        pass
    d=defaultdict(list)
    for lis in orbital_window_array:
        d[lis[0]].append(lis[1],)
    orbital_window_array=[list(x for y in i for x in y) for i in d.items()]                
    tddft_keywords_array=sp_keywords_array
    for m in range(0,len(orbital_window_array)):
        if orbital_window_array[m][0]=='N':
            tddft_input_file=path_out_LES+r"\TDDFT_%s-edge.inp" %orbital_window_array[m][0]
            tddft_output_file=path_out_LES+r"\TDDFT_%s-edge.out" %orbital_window_array[m][0]    
            tddft_error_file=path_out_LES+r"\TDDFT_%s-edge_error.txt" %orbital_window_array[m][0]
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
                for item in geom_array:
                    tddft_input.writelines(["%s " %item])
                tddft_input.close()
            #run tddft calculations
            tddft_out=open(tddft_output_file, "w") 
            tddft_err=open(tddft_error_file,"w") 
            p=sp.Popen([ORCA,tddft_input_file], stdout=tddft_out, stderr=tddft_err)
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
    
    sp_time=(sp_timer-start)/60
    tddft_time=(stop-sp_timer)/60
    
    print("\tSingle Point time is: "+ str(round(sp_time,3)) + " minutes")
    print("\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")
    
    log_file.write("\n\tSingle Point Energy time is: "+ str(round(sp_time,3)) + " minutes")
    log_file.write("\n\tTime-Dependent Density Functional Theory calculation time is: "+ str(round(tddft_time,3)) + " minutes")
    
    log_file.close()

    #run LC1
    p=sp.Popen(['python','U:\XAS_analysis\LC1\LC1.py',
        'U:\XAS_analysis\E2\Imdz_Solution_35C_APS.txt','0','3','7','-offset','38',
        str(path_out_LES+'\TDDFT_N-edge.out'), 
        'U:\XAS_analysis\E2\E2_Imdz_Solution_35C_APS_2018-08-20_16-33-03\Imdz_Solution_35C_APS.txt_fitted_peaks_param.txt'],
        stdout=sp.PIPE, stderr=sp.PIPE)
    p_status=p.wait()
    output, err = p.communicate()
    Norm_Trans_theory_data=path_in_LES+r'\..\LC1\NormTranslatedTheoryData.txt'
    theory_xdata=[]
    theory_ydata=[] 
    with open (Norm_Trans_theory_data,'r') as theory_file:
        lines=theory_file.readlines()
        for x in lines:
            theory_xdata.append(x.split('\t')[0])
            theory_ydata.append((x.split('\t')[1]).replace("\n",""))
    theory_file.close() 
    Exp_data=path_in_LES+r'\..\E2\Imdz_Solution_35C_APS.txt'
    exp_xdata=[]
    exp_ydata=[] 
    with open (Exp_data,'r') as exp_file:
        lines=exp_file.readlines()[38:]
        for line in lines:
            exp_xdata.append(line.split('  ')[1])
            exp_ydata.append(line.split('  ')[5])
           
