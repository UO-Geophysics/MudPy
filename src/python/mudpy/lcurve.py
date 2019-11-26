#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:19:39 2019

@author: degoldbe
"""

# Create L-curve to select regularization parameter
# take in home, project directory, run_name, scan .log files for Lm and L2 and create L curve
from glob import glob
import matplotlib.pyplot as plt
from numpy import zeros,c_

home='/Users/degoldbe/SlipInversion/'
project_name='OLD_2019RIDGECREST/Ridgecrest2019-Finale'
run_names=['2xM6_M7_SM_GPS_INSAR_T71Dx10_T64Ax6','FINE_2xM6_M7_SM_GPS_INSAR_T71Dx10_T64Ax6']
log_files=[]

for r in range(len(run_names)):
    logs=glob(home+project_name+'/output/inverse_models/models/'+run_names[r]+'.????.log')
    if r==0:
        log_files=logs
    else:
        log_files=log_files+logs

N=len(log_files)

sol_norm=zeros(N)
res_norm=zeros(N)
regularization=zeros(N)
M0s=zeros(N)

for r in range(len(log_files)):
    f=open(log_files[r],'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'lambda_spatial' in line:
            lambda_spatial=line.split('=')[-1]
        if 'L2 =' in line:
            L2=line.split('=')[-1]
        if 'Lm =' in line:
            Lm=line.split('=')[-1]
        if 'M0(N-m) =' in line:
            M0=line.split('=')[-1]
            loop_go=False
            f.close()
            print M0
    sol_norm[r]=Lm
    res_norm[r]=L2
    lambda_spatial="{0:.2f}".format(float(lambda_spatial))
    M0="{0:.2E}".format(float(M0))
    regularization[r]=lambda_spatial
    M0s[r]=M0


fig, ax = plt.subplots()
ax1=plt.subplot(111)
ax1.scatter(res_norm,sol_norm)
for i, txt in enumerate(regularization):
    ax1.annotate('$\lambda$='+str(txt)+'; M$_0$='+str(M0s[i]), (res_norm[i]+.002, sol_norm[i]+.002),fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1,100])
plt.xlim([.03,.3])
plt.xlabel('Residual Norm ||Gm$_\lambda$-d||$_2$')
plt.ylabel('Solution Norm ||m$_\lambda$||$_2$')
#ax2=plt.subplot(122)
#ax2.scatter(res_norm,sol_norm)
#for i,txt in enumerate(M0s):
#    ax2.annotate('M0='+str(txt), (res_norm[i]+.002, sol_norm[i]+.002),fontsize=12)
#    #ax.annotate(txt, (res_norm[i], sol_norm[i]))
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([1,100])
#plt.xlim([.03,.3])
#plt.xlabel('Residual Norm ||Gm$_\lambda$-d||$_2$')
#plt.ylabel('Solution Norm ||m$_\lambda$||$_2$')

f=open('/Users/degoldbe/Desktop/Lcurve.txt','w')
for r in range(len(res_norm)):
    f.write(str(sol_norm[r])+'\t'+str(res_norm[r])+'\t'+str(regularization[r])+'\t'+str(M0s[r])+'\n')
f.close()