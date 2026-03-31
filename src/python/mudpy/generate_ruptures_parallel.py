#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 10:45:24 2019

@author: degoldbe
"""

def run_parallel_generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
        load_distances,distances_name,UTM_zone,tMw,model_name,hurst,Ldip,Lstrike,
        num_modes,Nrealizations,rake,rise_time,rise_time_depths0,rise_time_depths1,time_epi,max_slip,
        source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,force_magnitude,
        force_area,mean_slip_name,hypocenter,slip_tol,force_hypocenter,
        no_random,use_hypo_fraction,shear_wave_fraction_shallow,shear_wave_fraction_deep,
        max_slip_rule,rank,size):
    
    '''
    Depending on user selected flags parse the work out to different functions
    '''
    
    from numpy import load,save,genfromtxt,log10,cos,sin,deg2rad,savetxt,zeros,where,argmin
    from time import gmtime, strftime
    from numpy.random import shuffle
    from mudpy import fakequakes
    from obspy import UTCDateTime
    from obspy.taup import TauPyModel
    import geopy.distance
    import warnings

    #I don't condone it but this cleans up the warnings
    warnings.filterwarnings("ignore")

    # Fix input formats
    rank=int(rank)
    size=int(size)
    if time_epi=='None':
        time_epi=None
    else:
        time_epi=UTCDateTime(time_epi)
    rise_time_depths=[rise_time_depths0,rise_time_depths1]
    tMw=tMw.split(',')
    target_Mw=zeros(len(tMw))
    for rMw in range(len(tMw)):
        target_Mw[rMw]=float(tMw[rMw])
        

    #Should I calculate or load the distances?
    if load_distances==1:  
        Dstrike=load(home+project_name+'/data/distances/'+distances_name+'.strike.npy')
        Ddip=load(home+project_name+'/data/distances/'+distances_name+'.dip.npy')
    else:
        Dstrike,Ddip=fakequakes.subfault_distances_3D(home,project_name,fault_name,slab_name,UTM_zone)
        save(home+project_name+'/data/distances/'+distances_name+'.strike.npy',Dstrike)
        save(home+project_name+'/data/distances/'+distances_name+'.dip.npy',Ddip)
    

    #Read fault and prepare output variable
    whole_fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    
    #Get structure model
    vel_mod_file=home+project_name+'/structure/'+model_name
    
    #Get TauPyModel
    velmod = TauPyModel(model=home+project_name+'/structure/'+model_name.split('.')[0])
    
    # Define the subfault hypocenter (if hypocenter is prescribed)
    if hypocenter is None:
        shypo=None
    else:
        dist=((whole_fault[:,1]-hypocenter[0])**2+(whole_fault[:,2]-hypocenter[1])**2)**0.5
        shypo=argmin(dist)
        
        # Re-define hypocenter as the coordinates of the hypocenter subfault in
        # case the original hypocenter did not perfectly align with a subfault
        hypocenter = whole_fault[shypo,1:4]
    

    #Now loop over the number of realizations
    realization=0
    if rank==0:
        print('Generating rupture scenarios')
    for kmag in range(len(target_Mw)):
        if rank==0:
            print('... Calculating ruptures for target magnitude Mw = '+str(target_Mw[kmag]))
        for kfault in range(Nrealizations):
            if kfault%1==0 and rank==0:
                print('... ... working on ruptures '+str(ncpus*realization)+' to ' + str(ncpus*(realization+1)-1) + ' of '+str(Nrealizations*size*len(target_Mw)))
                #print '... ... working on ruptures '+str(ncpus*realization+rank)+' of '+str(Nrealizations*size-1)
            
            #Prepare output
            fault_out=zeros((len(whole_fault),15))
            fault_out[:,0:8]=whole_fault[:,0:8]
            fault_out[:,10:12]=whole_fault[:,8:]   
            
            #Sucess criterion
            success=False
            while success==False:
                #Select only a subset of the faults based on magnitude scaling
                current_target_Mw=target_Mw[kmag]
                ifaults,hypo_fault,Lmax,Wmax,Leff,Weff,option,Lmean,Wmean=fakequakes.select_faults(whole_fault,Dstrike,Ddip,current_target_Mw,num_modes,scaling_law,
                                    force_area,no_shallow_epi=False,no_random=no_random,subfault_hypocenter=shypo,use_hypo_fraction=use_hypo_fraction)
                fault_array=whole_fault[ifaults,:]
                Dstrike_selected=Dstrike[ifaults,:][:,ifaults]
                Ddip_selected=Ddip[ifaults,:][:,ifaults]
                
                #Determine correlation lengths from effective length.width Leff and Weff
                if Lstrike=='MB2002': #Use scaling
                    #Ls=10**(-2.43+0.49*target_Mw)
                    Ls=2.0+(1./3)*Leff
                elif Lstrike=='auto':
                    Ls=17.7+0.34*Leff
                else:
                    Ls=Lstrike
                if Ldip=='MB2002': #Use scaling
                    #Ld=10**(-1.79+0.38*target_Mw)
                    Ld=1.0+(1./3)*Weff
                elif Ldip=='auto':
                    Ld=6.8+0.4*Weff
                else:
                    Ld=Ldip
                
                #Get the mean uniform slip for the target magnitude
                if mean_slip_name==None:
                    mean_slip,mu=fakequakes.get_mean_slip(target_Mw[kmag],fault_array,vel_mod_file)
                else:
                    foo,mu=fakequakes.get_mean_slip(target_Mw[kmag],fault_array,vel_mod_file)
                    mean_fault=genfromtxt(mean_slip_name)
                    mean_slip=(mean_fault[:,8]**2+mean_fault[:,9]**2)**0.5
                    
                    #keep onlt faults that have man slip inside the fault_array seelcted faults
                    mean_slip=mean_slip[ifaults]
                    
                    #get the area in those selected faults
                    area=fault_array[:,-2]*fault_array[:,-1]
                    
                    #get the moment in those selected faults
                    moment_on_selected=(area*mu*mean_slip).sum()
                    
                    #target moment
                    target_moment=10**(1.5*target_Mw[kmag]+9.1)
                    
                    #How much do I need to upscale?
                    scale_factor=target_moment/moment_on_selected
                    
                    #rescale the slip
                    mean_slip = mean_slip*scale_factor
                    
                    
                    #Make sure mean_slip has no zero slip faults
                    izero=where(mean_slip==0)[0]
                    mean_slip[izero]=slip_tol
                
                #Get correlation matrix
                C=fakequakes.vonKarman_correlation(Dstrike_selected,Ddip_selected,Ls,Ld,hurst)
                
                # Lognormal or not?
                if lognormal==False:
                    #Get covariance matrix
                    C_nonlog=fakequakes.get_covariance(mean_slip,C,target_Mw[kmag],fault_array,vel_mod_file,slip_standard_deviation) 
                    #Get eigen values and eigenvectors
                    eigenvals,V=fakequakes.get_eigen(C_nonlog)
                    #Generate fake slip pattern
                    rejected=True
                    while rejected==True:
#                        slip_unrectified,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=kfault)
                        slip_unrectified,success=fakequakes.make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=None)
                        slip,rejected,percent_negative=fakequakes.rectify_slip(slip_unrectified,percent_reject=13)
                        if rejected==True:
                            print('... ... ... negative slip threshold exceeeded with %d%% negative slip. Recomputing...' % (percent_negative))
                else:
                    #Get lognormal values
                    C_log,mean_slip_log=fakequakes.get_lognormal(mean_slip,C,target_Mw[kmag],fault_array,vel_mod_file,slip_standard_deviation)               
                    #Get eigen values and eigenvectors
                    eigenvals,V=fakequakes.get_eigen(C_log)
                    #Generate fake slip pattern
#                    slip,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=kfault)
                    slip,success=fakequakes.make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=None)
            
                #Slip pattern sucessfully made, moving on.
                #Rigidities
                foo,mu=fakequakes.get_mean_slip(target_Mw[kmag],whole_fault,vel_mod_file)
                fault_out[:,13]=mu
                
                #Calculate moment and magnitude of fake slip pattern
                M0=sum(slip*fault_out[ifaults,10]*fault_out[ifaults,11]*mu[ifaults])
                Mw=(2./3)*(log10(M0)-9.1)
                
                #Check max_slip_rule
                if max_slip_rule==True:
                    
                    max_slip_from_rule=10**(-4.94+0.71*Mw) #From Allen & Hayes, 2017
                    max_slip_tolerance = 3
                    
                    if slip.max() > max_slip_tolerance*max_slip_from_rule:
                        success = False
                        print('... ... ... max slip condition violated max_slip_rule, recalculating...')
                
                #Force to target magnitude
                if force_magnitude==True:
                    M0_target=10**(1.5*target_Mw[kmag]+9.1)
                    M0_ratio=M0_target/M0
                    #Multiply slip by ratio
                    slip=slip*M0_ratio
                    #Recalculate
                    M0=sum(slip*fault_out[ifaults,10]*fault_out[ifaults,11]*mu[ifaults])
                    Mw=(2./3)*(log10(M0)-9.1)
                    
                #check max_slip again
                if slip.max() > max_slip:
                    success=False
                    print('... ... ... max slip condition violated due to force_magnitude=True, recalculating...')
            
            
            #Get stochastic rake vector
            stoc_rake=fakequakes.get_stochastic_rake(rake,len(slip))
            
            #Place slip values in output variable
            fault_out[ifaults,8]=slip*cos(deg2rad(stoc_rake))
            fault_out[ifaults,9]=slip*sin(deg2rad(stoc_rake))
            
            #Move hypocenter to somewhere with a susbtantial fraction of peak slip
#            slip_fraction=0.25
#            islip=where(slip>slip.max()*slip_fraction)[0]
#            shuffle(islip) #randomize
#            hypo_fault=ifaults[islip[0]] #select first from randomized vector
            
            #Calculate and scale rise times
            rise_times=fakequakes.get_rise_times(M0,slip,fault_array,rise_time_depths,stoc_rake,rise_time,option=option)
            
            #Place rise_times in output variable
            fault_out[:,7]=0
            fault_out[ifaults,7]=rise_times
            
            #Calculate rupture onset times
            if force_hypocenter==False: #Use random hypo, otehrwise force hypo to user specified
                hypocenter=whole_fault[hypo_fault,1:4]

            
            # edit ...
            # if rise_time==2:
            #     shear_wave_fraction_shallow=1/60/60/24*2
            #     shear_wave_fraction_deep=   1/60/60/24*2
            #     print("c")
            # else:
            #     L_rescale=(Lmax/1500)+2.5
            #     shear_wave_fraction_shallow=1/60/60/24*2*(3.5/L_rescale)**2
            #     shear_wave_fraction_deep=   1/60/60/24*2*(3.5/L_rescale)**2
            #     print("m")
            if rise_time=='SSE':  #For moedling SSes
                shear_wave_fraction_shallow=1/60/60/24*2
                shear_wave_fraction_deep=   1/60/60/24*2
            else: #regular EQs, do nothing
                pass
            
            t_onset,length2fault=fakequakes.get_rupture_onset(home,project_name,slip,fault_array,model_name,hypocenter,shypo,rise_time_depths,
                                                 M0,velmod,shear_wave_fraction_shallow=shear_wave_fraction_shallow,
                                                 shear_wave_fraction_deep=shear_wave_fraction_deep)
            fault_out[:,12]=0
            fault_out[ifaults,12]=t_onset
            
            fault_out[:,14]=0
            fault_out[ifaults,14]=length2fault/t_onset
            
            
            #Calculate location of moment centroid
            centroid_lon,centroid_lat,centroid_z=fakequakes.get_centroid(fault_out)
            
            #Calculate average risetime
            rise = fault_out[:,7]
            avg_rise = np.mean(rise[np.where(rise>0)[0]])
            
            # Calculate average rupture velocity
            lon_array = fault_out[:,1]
            lat_array = fault_out[:,2]
            vrupt = []
            
            for i in range(len(fault_array)):
                if t_onset[i] > 0:
                    # r = geopy.distance.geodesic((hypocenter[1], hypocenter[0]), (lat_array[i], lon_array[i])).km
                    # vrupt.append(r/t_onset[i])
                    vrupt.append(length2fault[i]/t_onset[i])
            
            avg_vrupt = np.mean(vrupt)
            
            #Write to file
            run_number=str(ncpus*realization+rank).rjust(6,'0')
            outfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.rupt'
            #                             'No,   lon,    lat,  z(km),strike,   dip, rise,  dura,  ss(m),ds(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)
            savetxt(outfile,fault_out,fmt='%d\t%10.6f\t%10.6f\t%8.4f\t%7.2f\t%7.2f\t%4.1f\t%.9e\t%.4e\t%.4e\t%10.2f\t%10.2f\t%.9e\t%.6e\t%.6e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa),velocity(km/s)')
            # savetxt(outfile,fault_out,fmt='%d\t%10.6f\t%10.6f\t%8.4f\t%7.2f\t%7.2f\t%4.1f\t%5.2f\t%5.2f\t%5.2f\t%10.2f\t%10.2f\t%5.2f\t%.6e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')

            #Write log file
            logfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
            f=open(logfile,'w')
            f.write('Scenario calculated at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+' GMT\n')
            f.write('Project name: '+project_name+'\n')
            f.write('Run name: '+run_name+'\n')
            f.write('Run number: '+run_number+'\n')
            f.write('Velocity model: '+model_name+'\n')
            f.write('No. of KL modes: '+str(num_modes)+'\n')
            f.write('Hurst exponent: '+str(hurst)+'\n')
            f.write('Corr. length used Lstrike: %.2f km\n' % Ls)
            f.write('Corr. length used Ldip: %.2f km\n' % Ld)
            f.write('Slip std. dev.: %.3f km\n' % slip_standard_deviation)
            f.write('Maximum length Lmax: %.2f km\n' % Lmax)
            f.write('Maximum width Wmax: %.2f km\n' % Wmax)
            f.write('Effective length Leff: %.2f km\n' % Leff)
            f.write('Effective width Weff: %.2f km\n' % Weff)
            f.write('Target magnitude: Mw %.4f\n' % target_Mw[kmag])
            f.write('Actual magnitude: Mw %.4f\n' % Mw)
            f.write('Hypocenter (lon,lat,z[km]): (%.6f,%.6f,%.2f)\n' %(hypocenter[0],hypocenter[1],hypocenter[2]))
            f.write('Hypocenter time: %s\n' % time_epi)
            f.write('Centroid (lon,lat,z[km]): (%.6f,%.6f,%.2f)\n' %(centroid_lon,centroid_lat,centroid_z))
            f.write('Source time function type: %s\n' % source_time_function)
            f.write('Average Risetime (s): %.9e\n' % avg_rise)
            f.write('Average Rupture Velocity (km/s): %.9e\n' % avg_vrupt)
            f.write('Avg. length: %.2f km\n' % Lmean)
            f.write('Avg. width: %.2f km\n' % Wmean)
            f.write('Class: %d type' % option)
            # f.write('Average Risetime (s): %.2f\n' % avg_rise)
            # f.write('Average Rupture Velocity (km/s): %.2f' % avg_vrupt)
            f.close()
            
            realization+=1
            
            print('Run number: '+run_number+'\n')

#If main entry point
if __name__ == '__main__':
    import sys
    import numpy as np
    from mpi4py import MPI
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #Map command line arguments to funciton arguments:
    if sys.argv[1]=='run_parallel_generate_ruptures':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        run_name=sys.argv[4]
        fault_name=sys.argv[5]
        slab_name=sys.argv[6]
        if slab_name=='None':
            slab_name=None
        mesh_name=sys.argv[7]
        if mesh_name=='None':
            mesh_name=None
        load_distances=int(sys.argv[8])
        distances_name=sys.argv[9]
        UTM_zone=sys.argv[10]
        tMw=sys.argv[11]
        model_name=sys.argv[12]
        hurst=float(sys.argv[13])
        Ldip=sys.argv[14]
        Lstrike=sys.argv[15]
        num_modes=int(sys.argv[16])
        Nrealizations=int(sys.argv[17])
        rake=float(sys.argv[18])
        rise_time=sys.argv[19]
        rise_time_depths0=int(sys.argv[20])
        rise_time_depths1=int(sys.argv[21])
        time_epi=sys.argv[22]
        max_slip=float(sys.argv[23])
        source_time_function=sys.argv[24]
        lognormal=sys.argv[25]
        if lognormal=='True':
            lognormal=True
        elif lognormal=='False':
            lognormal=False
        slip_standard_deviation=float(sys.argv[26])
        scaling_law=sys.argv[27]        
        ncpus=int(sys.argv[28])
        force_magnitude=sys.argv[29]
        if force_magnitude=='True':
            force_magnitude=True
        elif force_magnitude=='False':
            force_magnitude=False
        force_area=sys.argv[30]
        if force_area=='True':
            force_area=True
        elif force_area=='False':
            force_area=False
        mean_slip_name=sys.argv[31]
        if mean_slip_name == 'None':
            mean_slip_name=None
        hypocenter=sys.argv[32]
        if hypocenter == 'None':
            hypocenter=None
        else:
            hypocenter_lon=float(hypocenter.split(',')[0].split('[')[1])
            hypocenter_lat=float(hypocenter.split(',')[1])
            hypocenter_dep=float(hypocenter.split(',')[2].split(']')[0])
            hypocenter=np.array([hypocenter_lon,hypocenter_lat,hypocenter_dep])
        slip_tol=float(sys.argv[33])
        force_hypocenter=sys.argv[34]
        if force_hypocenter=='True':
            force_hypocenter=True
        elif force_hypocenter=='False':
            force_hypocenter=False
        no_random=sys.argv[35]
        if no_random=='True':
            no_random=True
        elif no_random=='False':
            no_random=False
        use_hypo_fraction=sys.argv[36]
        if use_hypo_fraction=='True':
            use_hypo_fraction=True
        if use_hypo_fraction=='False':
            use_hypo_fraction=False
        shear_wave_fraction_shallow=float(sys.argv[37])
        shear_wave_fraction_deep=float(sys.argv[38])
        max_slip_rule=sys.argv[39]
        if max_slip_rule=='True':
            max_slip_rule=True
        if max_slip_rule=='False':
            max_slip_rule=False
        
        run_parallel_generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
        load_distances,distances_name,UTM_zone,tMw,model_name,hurst,Ldip,Lstrike,
        num_modes,Nrealizations,rake,rise_time,rise_time_depths0,rise_time_depths1,time_epi,max_slip,
        source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,force_magnitude,
        force_area,mean_slip_name,hypocenter,slip_tol,force_hypocenter,
        no_random,use_hypo_fraction,shear_wave_fraction_shallow,shear_wave_fraction_deep,
        max_slip_rule,rank,size)
    else:
        print("ERROR: You're not allowed to run "+sys.argv[1]+" from the shell or it does not exist")
        
