#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:22:47 2022

@author: dmelgarm
"""

from numpy import genfromtxt
from os import environ,path
from obspy.taup import taup_create, TauPyModel
import os

#mudpy source folder
vel_mod_file = '/Users/dmelgarm/FakeQuakes/ONC_debug/structure/cascadia2.mod'
home='/Users/dmelgarm/FakeQuakes/'
project_name='ONC_debug'
background_model='PREM'
#load user specified .mod infromation
structure = genfromtxt(vel_mod_file)

#load background velocity structure
if background_model=='PREM':
    
    bg_model_file='/Users/dmelgarm/code/MudPy/src/aux/prem.nd'
    
    #Q values
    Qkappa=1300
    Qmu=600
    
    #Write new _nd file one line at a time
    nd_name=path.basename(vel_mod_file).split('.')[0]
    nd_name=nd_name+'.nd'
    f=open(home+project_name+'/structure/'+nd_name,'w')
    
    #initalize
    ztop=0
    
    for k in range(len(structure)-1):
        
        #Variables for writing to file
        zbot=ztop+structure[k,0]
        vp=structure[k,2]
        vs=structure[k,1]
        rho=structure[k,3]
        
        # Write to the file
        line1=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (ztop,vp,vs,rho,Qkappa,Qmu))
        line2=('%8.2f\t%8.5f   %7.5f   %7.5f\t%6.1f     %5.1f\n' % (zbot,vp,vs,rho,Qkappa,Qmu))
        f.write(line1)
        f.write(line2)
        
        #update
        ztop=zbot

    
    #now read PREM file libe by libne and find appropriate depth tos tart isnerting
    fprem=open(bg_model_file,'r')
    found_depth=False
    
    while True:
        
        line=fprem.readline()
        
        if line=='': #End of file
            break
        
        if found_depth==False:
            #Check that it's not a keyword line like 'mantle'
            if len(line.split())>1:
                
                #not a keyword, waht depth are we at?
                current_depth=float(line.split()[0])
                
                if current_depth > zbot: #PREM depth alrger than .mod file last line
                    found_depth=True
                    f.write('mantle\n')
                
        #Ok you have found the depth write it to file
        if found_depth == True:
            f.write(line)
            
    
    fprem.close()
    f.close()

    # make TauPy npz
    taup_in=home+project_name+'/structure/'+nd_name
    taup_out=home+project_name+'/structure/'
    
    # taup_create.build_taup_model(taup_in,output_folder=taup_out)
    
    model_name = os.path.splitext(os.path.basename(taup_in))[0]
    output_filename = os.path.join(taup_out, model_name + ".npz")

    mod_create = taup_create.TauPCreate(input_filename=taup_in,
                            output_filename=output_filename)
    
    mod_create.min_delta_p = 0.001
    mod_create.max_delta_p = 0.1
    mod_create.load_velocity_model()
    mod_create.run()
    
    
    
    velmod = TauPyModel(model='/Users/dmelgarm/FakeQuakes/ONC_debug/structure/cascadia2.npz')
    
    
    
    
    
    #testing stuff
    zs=13.866199999999999
    dist_in_degs=6.843948273047394
    Spaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['S','s'])
    Spaths.plot_rays(plot_type='cartesian')