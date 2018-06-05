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
import datetime
import os
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
from collections import OrderedDict
import pandas as pd

filename=r"N1s_Imidazole_ISEELS_Hitchcock_Norm_Athena.txt" #provided by user 

#resultsdir = r"E2_"+filename+r"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#working_dir=os.getcwd()
working_dir=r"C:\Users\dmg81179\Desktop\Code_Development\2018-June_peak_fitting_E2\working_dir"
path_in=working_dir
path_out=working_dir #+r"\\"+resultsdir
#os.makedirs(path_out)
#path_out=working_dir

edge_data_table=path_in+r"\edge_data.txt"
data_file=path_in+r"\%s" %filename 
fitted_peak=path_out+r"\fitted_peaks"
peak_params=path_out+r"\peak_params.txt"


edge_data=np.array([])
data=np.array([])


with open(edge_data_table,"r") as edge_data_table:
    next (edge_data_table)
    next (edge_data_table)
    for line in edge_data_table:
        line=line.split()
        edge_data=np.append(edge_data,line)
    edge_data=edge_data.reshape(int(len(edge_data)/5),5)
edge_data_table.close()

with open (data_file) as data_file:
    lines_after_heading=data_file.readlines()[38:]# 38 is default but it can change 
    for line in lines_after_heading:
        line=line.split()
        data=np.append(data,line)
    data=data.reshape(int(len(data)/6),6)#6 is the number of columns and it can change 
data_file.close()

#xdata and ydata
ycol=2
xdata=data[:,0].astype(float) #xdata is always the first column
ydata=data[:,ycol].astype(float) # ask user for the column where the ydata is 

#determine e0 
dif1=np.diff(ydata)/np.diff(xdata)
e0_pos=np.where(dif1==max(dif1))[0][0]
e0=xdata[e0_pos]
e_edge=abs(edge_data[:,4].astype(float)-e0)
e1_pos=np.where(e_edge==min(e_edge))[0][0]
element=edge_data[e1_pos][0]

"""
#This was copied from the matlab legacy code and was not used
#determining ranges of pre-edge line and post-edge poly.
fit_limits=[150,60,0,0] #input from user
pre_pos_i_ls=abs(xdata-e0+fit_limits[0])
pre_pos_i=np.where(pre_pos_i_ls==min(pre_pos_i_ls))[0][0]
pre_pos_f_ls=abs(xdata-e0+fit_limits[1])
pre_pos_f=np.where(pre_pos_f_ls==min(pre_pos_f_ls))[0][0]

post_pos_i_ls=abs(xdata-e0+fit_limits[2])
post_pos_i=np.where(post_pos_i_ls==min(post_pos_i_ls))[0][0]
post_pos_f=np.size(xdata)-fit_limits[3]-1
"""
#determine fitting range
fit_range=[5,5]
fit_pos_i_ls=abs(xdata-e0+fit_range[0])
fit_pos_i=np.where(fit_pos_i_ls==min(fit_pos_i_ls))[0][0]
fit_pos_f_ls=abs(xdata-e0-fit_range[1])
fit_pos_f=np.where(fit_pos_f_ls==min(fit_pos_f_ls))[0][0]
fit_xdata=xdata[fit_pos_i:fit_pos_f].astype(float)
fit_ydata=ydata[fit_pos_i:fit_pos_f].astype(float)

#functions used for fitting
#from lmfit.models import StepModel, PseudoVoigtModel, GaussianModel
from lmfit.models import GaussianModel, StepModel

#determination of number of gaussian and their initial position guesses
smooth_ydata=scipy.signal.savgol_filter(fit_ydata,11,2)
dif2=np.diff(smooth_ydata)/np.diff(fit_xdata)
e0_pos2=np.where(dif2==max(dif2))[0][0]
posit=np.array((np.where((dif2[1:]<0)*(dif2[0:-1]>0))),dtype='int')+1
posit=np.unique(np.sort(np.append(posit,np.array((np.where((dif2[1:]>0)*(dif2[0:-1]<0))),dtype='int')+1)))
posit=posit[posit>e0_pos2]
funccenter=fit_xdata[posit]


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


##First fitting attempt 
gaussnum=len(funccenter)
funcnum=gaussnum
#initial guess x0F, lower bound lb, upper bound up
x0f=np.zeros((funcnum,4))
lb=np.zeros((funcnum,4))
ub=np.zeros((funcnum,4))
#initial guess for error function
step1=StepModel(form='arctan', prefix='step1_')
#pars.update(step2.guess(y,x=x))
pars=step1.make_params()
pars['step1_center'].set(e0+4, min=e0+3, max=e0+6)
pars['step1_amplitude'].set(0.5, min=0.1, max=1)
pars['step1_sigma'].set(0.5, min=0.3, max=0.8)
mod=step1
#initial guess for gaussian
x0f [:funcnum,2]=funccenter[:]
for n0 in range (0,funcnum):
    x0f[n0,0:2]=[0.5,0.5]
    lb[n0,:]=[0.2,0.2,x0f[n0,2]-1,0]
    ub[n0,:]=[3,0.7,x0f[n0,2]+1,0.1]
    gauss=GaussianModel(prefix='g%s_'%int(n0+1))
    pars.update(gauss.make_params())
    pars['g%s_amplitude'%int(n0+1)].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
    pars['g%s_sigma'%int(n0+1)].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
    pars['g%s_center'%int(n0+1)].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
    mod+=gauss          
#fitting the functions
init = mod.eval(pars, x=fit_xdata)
out = mod.fit(fit_ydata, pars, x=fit_xdata)
Y_fitted=out.best_fit
R_sqr=1 - out.residual.var() / np.var(fit_ydata)

#reading fitted peaks parameters from first fitting attempt
v=[]
for param in out.params.values():
    v.append("%s:  %f" % (param.name, param.value))
    #param_values.append("%s:  %f +/- %f (init = %f)" % (param.name, param.value, param.stderr, param.init_value))
param_values=[]
for item in v:
    param_values.append(item.split(":"))
function_keywords=["step"]
d=OrderedDict([("step",[])])
for i in range(1,funcnum+1):
    function_keywords.append("g%s"%i)
    d["g%s"%i]=[]
param_keywords=["amplitude","sigma","center","fwhm","height",]
for n in function_keywords:
    for c in range (0,len(param_values)):
        if n in param_values[c][0]:
            d[n].append(param_values[c][1])
        c+=1
funccenter_new=np.array([])
func_diff=[]
for n in range(1,gaussnum):
    func_diff.append(abs(float(d['g%s'%str(n+1)][2])-float(d['g%s'%str(n)][2])))
    if abs(float(d['g%s'%str(n+1)][2])-float(d['g%s'%str(n)][2])) >0.5:
        funccenter_new=np.append(funccenter_new,float(d['g%s'%str(n+1)][2]))
        funccenter_new=np.append(funccenter_new,float(d['g%s'%str(n)][2]))
        funccenter_new=np.unique(np.sort(funccenter_new))
if any(t<0.5 for t in func_diff):
    # second fitting attempt, after peaks with same energy position removed
    gaussnum=len(funccenter_new)
    funcnum=gaussnum
    #initial guess x0F, lower bound lb, upper bound up
    x0f=np.zeros((funcnum,4))
    lb=np.zeros((funcnum,4))
    ub=np.zeros((funcnum,4))
    #initial guess for error function
    step1=StepModel(form='arctan', prefix='step1_')
    #pars.update(step2.guess(y,x=x))
    pars=step1.make_params()
    pars['step1_center'].set(e0+4, min=e0+3, max=e0+6)
    pars['step1_amplitude'].set(0.5, min=0.1, max=1)
    pars['step1_sigma'].set(0.5, min=0.3, max=0.8)
    mod=step1
    #initial guess for gaussian
    x0f [:funcnum,2]=funccenter_new[:]
    for n0 in range (0,funcnum):
        x0f[n0,0:2]=[0.5,0.5]
        lb[n0,:]=[0.2,0.2,x0f[n0,2]-1,0]
        ub[n0,:]=[3,0.7,x0f[n0,2]+1,0.1]
        gauss=GaussianModel(prefix='g%s_'%int(n0+1))
        pars.update(gauss.make_params())
        pars['g%s_amplitude'%int(n0+1)].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
        pars['g%s_sigma'%int(n0+1)].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
        pars['g%s_center'%int(n0+1)].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
        mod+=gauss          
    #fitting the functions
    init = mod.eval(pars, x=fit_xdata)
    out = mod.fit(fit_ydata, pars, x=fit_xdata)
    Y_fitted=out.best_fit
    R_sqr=1 - out.residual.var() / np.var(fit_ydata)
    
    #reading fitted peaks parameters from second fitting attempt
    v=[]
    for param in out.params.values():
        v.append("%s:  %f" % (param.name, param.value))
        #param_values.append("%s:  %f +/- %f (init = %f)" % (param.name, param.value, param.stderr, param.init_value))
    param_values=[]
    for item in v:
        param_values.append(item.split(":"))
    function_keywords=["step"]
    d=OrderedDict([("step",[])])
    for i in range(1,funcnum+1):
        function_keywords.append("g%s"%i)
        d["g%s"%i]=[]
    param_keywords=["amplitude","sigma","center","fwhm","height",]
    for n in function_keywords:
        for c in range (0,len(param_values)):
            if n in param_values[c][0]:
                d[n].append(param_values[c][1])
            c+=1

#plot results
fig=plt.figure() 
plot_components = True
plt.plot(xdata, ydata, 'b',label='Experimental Data')
#plt.plot(fit_xdata, init, 'k--') #plot with initial guess
plt.plot(fit_xdata, out.best_fit, 'r-', label='Fitted data')

if plot_components:
    comps = out.eval_components(x=fit_xdata)
    plt.plot(fit_xdata,comps['step1_'], 'k--',label='Step Function')
    for x in range (1,gaussnum+1):
        plt.plot(fit_xdata, comps['g%s_'%x], 'k--',label='Gaussian%s Function'%x)
plt.legend(loc='upper left')
ax = plt.gca()
ax.set_xlim([fit_xdata[0],fit_xdata[-1]+1])
ax.set_ylim([0,max(fit_ydata)+0.5])
plt.xlabel('Energy/ eV')
plt.ylabel('Intensity')
fig.show()
fig.savefig('Fitted_peaks.tif')
print("Goodness of fit (R-sqaured) is: %s" %R_sqr)

#write fitted peaks to output file 
fitted_peaks=pd.concat([pd.DataFrame({'energy': fit_xdata}), pd.DataFrame(comps)], axis=1)
fitted_peaks.to_csv(path_out+r'\fitted_peaks.txt', index=False, sep='\t', header=True)

#Write fitted params to output file 
fitted_peaks_param=pd.concat([pd.DataFrame({'parameter': param_keywords}),pd.DataFrame(dict([(k,pd.Series(v)) for k,v in d.items()]))], axis=1)
fitted_peaks_param.to_csv(path_out+r'\fitted_peaks_param.txt',index=False, sep='\t', header=True)

stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 