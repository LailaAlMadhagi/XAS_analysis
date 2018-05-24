# -*- coding: utf-8 -*-
"""
Created on Tue May 22 17:16:41 2018

@author: Laila Al-Madhagi 
email: fy11lham@leeds.ac.uk 
"""

import timeit
start = timeit.default_timer()
import datetime
import os
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt


filename=r"1_5mM_kaucn4_h2o_1.nor" #provided by user 

resultsdir = r"E2_"+filename+r"_"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
working_dir=os.getcwd()
path_in=working_dir
path_out=working_dir+r"\\"+resultsdir
os.makedirs(path_out)


edge_data_table=path_in+r"\edge_data.txt"
edge_data=np.array([])
data_file=path_in+r"\%s" %filename 
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
    lines_after_heading=data_file.readlines()[31:]# 31 is default but it can change 
    for line in lines_after_heading:
        line=line.split()
        data=np.append(data,line)
    data=data.reshape(int(len(data)/7),7)#7 is the number of columns and it can change 
data_file.close()

#xdata and ydata
ycol=3
xdata=data[:,0].astype(float) #xdata is always the first column
ydata=data[:,ycol].astype(float) # ask user for the column where the ydata is 

#determine e0 
dif1=np.diff(ydata)/np.diff(xdata)
e0_pos=np.where(dif1==max(dif1))[0][0]
e0=xdata[e0_pos]
e_edge=abs(edge_data[:,4].astype(float)-e0)
e1_pos=np.where(e_edge==min(e_edge))[0][0]
element=edge_data[e1_pos][0]


#determining ranges of pre-edge line and post-edge poly.
fit_limits=[150,60,0,0] #input from user
pre_pos_i_ls=abs(xdata-e0+fit_limits[0])
pre_pos_i=np.where(pre_pos_i_ls==min(pre_pos_i_ls))[0][0]
pre_pos_f_ls=abs(xdata-e0+fit_limits[1])
pre_pos_f=np.where(pre_pos_f_ls==min(pre_pos_f_ls))[0][0]

post_pos_i_ls=abs(xdata-e0+fit_limits[2])
post_pos_i=np.where(post_pos_i_ls==min(post_pos_i_ls))[0][0]
post_pos_f=np.size(xdata)-fit_limits[3]-1

#determine fitting range
fit_range=[15,60]
fit_pos_i_ls=abs(xdata-e0+fit_range[0])
fit_pos_i=np.where(fit_pos_i_ls==min(fit_pos_i_ls))[0][0]
fit_pos_f_ls=abs(xdata-e0-fit_range[1])
fit_pos_f=np.where(fit_pos_f_ls==min(fit_pos_f_ls))[0][0]+1
fit_xdata=xdata[fit_pos_i:fit_pos_f].astype(float)
fit_ydata=ydata[fit_pos_i:fit_pos_f].astype(float)

#functions used for fitting
from lmfit.models import StepModel, PseudoVoigtModel, GaussianModel
"""
#Error Function:
def f0(x,xdata):
    return x(1)/2*(np.erf((xdata-x(3))*1.665/x(2))+1)
#pseudo Viogt function
def f1 (x,xdata):
    return x(4)*x(1)/3.14*(((x(2)/2))/((xdata-x(3))^2 +(x(2)/2)^2))+(1-x(4))*x(1)/(x(2)*(2*3.14)^.5)*np.exp(-0.5*((xdata-x(3))*2.355/x(2))^2)
#pseudo Viogt2 function
def f2 (x,xdata):
    return x(4)*x(1)/3.14*(((x(2)/2))/((xdata-x(3))^2 +(x(2)/2)^2))+(1-x(4))*x(1)/(x(2)*(2*3.14)^.5)*np.exp(-0.5*((xdata-x(3))*2.355/x(2))^2)
#Gaussian function
def f3 (x,xdata):
    return x(1)/(x(2)*(2*3.14)^.5)*np.exp(-0.5*((xdata-x(3))*2.355/x(2))^2)
"""

#determination of number of gaussian and their initial position guesses
smooth_ydata=scipy.signal.savgol_filter(fit_ydata,11,2)
dif2=np.diff(fit_ydata)/np.diff(fit_xdata)
e0_pos2=np.where(dif2==max(dif2))[0][0]
posit=np.array((np.where((dif2[1:]<0)*(dif2[0:-1]>0))),dtype='int')+1
posit=posit[posit>e0_pos2]
funccenter=fit_xdata[posit]

#remove extraneous peaks
a=np.diff(funccenter)
b=np.where(a>3.5)#the value 3.5 how to change 
funccenter=funccenter[b]

#Peak fitting
erfnum=1
pvoigtnum=1
pvoigtnum2=0
gaussnum=len(funccenter)-pvoigtnum2-pvoigtnum
funcnum=pvoigtnum+erfnum+pvoigtnum2+gaussnum

#initial guess x0F, lower bound lb, upper bound up
x0f=np.zeros((funcnum,4))
lb=np.zeros((funcnum,4))
ub=np.zeros((funcnum,4))


#initial guess for error function
for n0 in range(0,erfnum):
    if erfnum==1:
        x0f[n0,:]=[0.8,edge_data[e1_pos,3],e0+10,0]
        lb[n0,:]=[0.8,edge_data[e1_pos,3],e0,0]
        ub[n0,:]=[0.8001,10,e0+25,0.1]
        #step=locals()[step+str(n0)]        
        step=StepModel(form='erf', prefix='step%s_'%erfnum)
        pars=step.make_params()
        pars['step%s_amplitude'%erfnum].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
        pars['step%s_sigma'%erfnum].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
        pars['step%s_center'%erfnum].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
        mod=step
    elif erfnum==2 and n0==0:
        x0f[n0,:]=[0.3,5,e0+15,0]
        lb[n0,:]=[0,edge_data[e1_pos,3],e0,0]
        ub[n0,:]=[0.5,15,e0+60,0.1]
        step=StepModel(form='erf', prefix='step%s_'%str(erfnum-1))
        pars=step.make_params()
        pars['step%s_amplitude'%str(erfnum-1)].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
        pars['step%s_sigma'%str(erfnum-1)].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
        pars['step%s_center'%str(erfnum-1)].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
        mod=step
    elif erfnum==2 and n0==1:
        x0f[n0,:]=[0.1,0,e0+15,0]
        lb[n0,:]=[0,0,e0,0]
        ub[n0,:]=[0.2,0.1,e0+20,0.1]
        step=StepModel(form='erf', prefix='step%s_'%erfnum)
        pars=step.make_params()
        pars['step%s_amplitude'%erfnum].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
        pars['step%s_sigma'%erfnum].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
        pars['step%s_center'%erfnum].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
        mod+=step
#initial guess for pseudo Voigt
x0f [erfnum:erfnum+pvoigtnum,2]=funccenter[0]
x0f [erfnum:erfnum+pvoigtnum,3]=0.8
for n0 in range (erfnum,erfnum+pvoigtnum):
    x0f[n0,0:2]=[0.5,edge_data[e1_pos,3]]
    lb[n0,:]=[0,edge_data[e1_pos,3],x0f[n0,2]-2,0.9]
    ub[n0,:]=[15,20,x0f[n0,2]+2,1]
    pvoigt=PseudoVoigtModel(prefix='pv%s_'%pvoigtnum)
    pars.update(pvoigt.make_params())
    pars['pv%s_amplitude'%pvoigtnum].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
    pars['pv%s_sigma'%pvoigtnum].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
    pars['pv%s_center'%pvoigtnum].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
    mod+=pvoigt
  
#initial guess for pseudo Voigt2
x0f [erfnum+pvoigtnum:erfnum+pvoigtnum+pvoigtnum2,2]=funccenter[1]
for n0 in range (erfnum+pvoigtnum,erfnum+pvoigtnum+pvoigtnum2):
    x0f[n0,0:2]=[0.5,edge_data[e1_pos,3]]
    lb[n0,:]=[0,edge_data[e1_pos,3],x0f[n0,2],0.5]
    ub[n0,:]=[15,edge_data[e1_pos,3]+10,x0f[n0,2]+10,1]
    pvoigt=PseudoVoigtModel(prefix='pv%s_'%pvoigtnum2)
    pars.update(pvoigt.make_params())
    pars['pv%s_amplitude'%pvoigtnum2].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
    pars['pv%s_sigma'%pvoigtnum2].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
    pars['pv%s_center'%pvoigtnum2].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
    mod+=pvoigt
#initial guess for gaussian
x0f [erfnum+pvoigtnum+pvoigtnum2:funcnum,2]=funccenter[pvoigtnum+pvoigtnum2:]
for n0 in range (erfnum+pvoigtnum+pvoigtnum2,funcnum):
    x0f[n0,0:2]=[0.2,edge_data[e1_pos,3]]
    lb[n0,:]=[0,0.5,x0f[n0,2]-3,0]
    ub[n0,:]=[15,20,x0f[n0,2]+3,0.1]
    gauss=GaussianModel(prefix='g%s_'%int(n0-erfnum-pvoigtnum-pvoigtnum2+1))
    pars.update(gauss.make_params())
    pars['g%s_amplitude'%int(n0-erfnum-pvoigtnum-pvoigtnum2+1)].set(x0f[n0][0], min=lb[n0][0],max=ub[n0][0])
    pars['g%s_sigma'%int(n0-erfnum-pvoigtnum-pvoigtnum2+1)].set(x0f[n0][1], min=lb[n0][1],max=ub[n0][1])
    pars['g%s_center'%int(n0-erfnum-pvoigtnum-pvoigtnum2+1)].set(x0f[n0][2], min=lb[n0][2],max=ub[n0][2])
    mod+=gauss
    
#fitting final function
init = mod.eval(pars, x=fit_xdata)
out = mod.fit(fit_ydata, pars, x=fit_xdata)

plot_components = True

plt.plot(xdata, ydata, 'b')
#plt.plot(fit_xdata, init, 'k--') #plot with initial guess
plt.plot(fit_xdata, out.best_fit, 'r-')

if plot_components:
    comps = out.eval_components(x=fit_xdata)
    plt.plot(fit_xdata, comps['step1_'], 'k--')
    plt.plot(fit_xdata, comps['pv1_'], 'b--')
    plt.plot(fit_xdata, comps['g1_'], 'b--')
    plt.plot(fit_xdata, comps['g2_'], 'b--')
    plt.plot(fit_xdata, comps['g3_'], 'b--')
    plt.plot(fit_xdata, comps['g4_'], 'b--')
ax = plt.gca()
ax.set_xlim([11850,12000])
ax.set_ylim([0,1.6])
plt.show()

"""
#this block is weird, assign functions 
#Fchar=['x,xdata']
f00='+x0F[tmp,0]/2*(np.erf((xdata-x0F[tmp,2])*1.665/x0F[tmp,1])+1)'
f01='+x0F[tmp,3]*x0F[tmp,0]/3.14*((x0F[tmp,1]/2)/((xdata-x0F[tmp,2])**2 + (x0F[tmp,1]/2)**2))+(1-x0F[tmp,3])*x0F[tmp,0]/(x0F[tmp,1]*(2*3.14)**.5)*np.exp(-0.5*((xdata-x0F[tmp,2])*2.355/x0F[tmp,1])**2)'
f02='+x0F[tmp,3]*x0F[tmp,0]/3.14*((x0F[tmp,1]/2)/((xdata-x0F[tmp,2])**2 + (x0F[tmp,1]/2)**2))+(1-x0F[tmp,3])*x0F[tmp,1]/(x0F[tmp,1]*(2*3.14)**.5)*np.exp(-0.5*((xdata-x0F[tmp,2])*2.355/x0F[tmp,1])**2)'
f03='+x0F[tmp,0]/(x0F[tmp,1]*(2*3.14)**.5)*np.exp(-0.5*((xdata-x0F[tmp,2])*2.355/x0F[tmp,1])**2)'
for n in range(1,funcnum+1):
    if n <= erfnum:
        mod=step
        
        if n==1:
            Fchar=[a.replace('tmp',str(n-1)) for a in Fchar]
        if n==2:
            Fchar=[a.replace('x(tmp,1)','(0.8-x(0,1)+x(1,1))') for a in Fchar]
        
    elif n > erfnum and n<=erfnum+pvoigtnum:
        mod+=pvoigt1
    elif n > erfnum+pvoigtnum and n <= erfnum+pvoigtnum+pvoigtnum2:
        mod+=pvoigt2
    else:
        mod+=gauss



def make_func(Fchar):
    funcstr='''\
def F(x,xdata):
    return {e}
    '''.format(e=Fchar)
    exec(funcstr)
    return F
F=make_func(Fchar)

def errfunc(x,xdata,ydata):
    err=ydata-f(x,xdata)
    return err
"""


stop = timeit.default_timer()
running_time=(stop-start)/60
print ("Running time is: "+ str(round(running_time,3)) + "minutes") 