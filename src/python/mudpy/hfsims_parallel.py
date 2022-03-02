'''
Module for high frequency simulations in parallel
'''

def run_parallel_hfsims(home,project_name,rupture_name,N,M0,sta,sta_lon,sta_lat,component,model_name,
                        rise_time_depths0,rise_time_depths1,moho_depth_in_km,total_duration,
                        hf_dt,stress_parameter,kappa,Qexp,Pwave,Swave,high_stress_depth,rank,size): 
    '''
    Run stochastic HF sims
    
    stress parameter is in bars
    '''
    
    from numpy import genfromtxt,pi,logspace,log10,mean,where,exp,arange,zeros,argmin,rad2deg,arctan2,real,savetxt,c_
    from pyproj import Geod
    from obspy.geodetics import kilometer2degrees
    from obspy.taup import TauPyModel
    from mudpy.forward import get_mu, read_fakequakes_hypo_time
    from mudpy import hfsims
    from obspy import Stream,Trace
    from sys import stdout
    from os import path,makedirs
    from mudpy.hfsims import is_subfault_in_smga
    import warnings

    rank=int(rank)
    
    if rank==0 and component=='N':
        #print out what's going on:
        out='''Running with input parameters:
        home = %s
        project_name = %s
        rupture_name = %s
        N = %s
        M0 (N-m) = %s
        sta = %s
        sta_lon = %s
        sta_lat = %s
        model_name = %s
        rise_time_depths = %s
        moho_depth_in_km = %s
        total_duration = %s
        hf_dt = %s
        stress_parameter = %s
        kappa = %s
        Qexp = %s
        component = %s
        Pwave = %s
        Swave = %s
        high_stress_depth = %s
        '''%(home,project_name,rupture_name,str(N),str(M0/1e7),sta,str(sta_lon),str(sta_lat),model_name,str([rise_time_depths0,rise_time_depths1]),
        str(moho_depth_in_km),str(total_duration),str(hf_dt),str(stress_parameter),
        str(kappa),str(Qexp),str(component),str(Pwave),str(Swave),str(high_stress_depth))
        print(out)

    if rank==0:
        out='''
        Rupture_Name = %s
        Station = %s
        Component (N,E,Z) = %s
        Sample rate = %sHz
        Duration = %ss
        '''%(rupture_name,sta,component,str(1/hf_dt),str(total_duration))
        print(out)
        
    #print 'stress is '+str(stress_parameter)

    #I don't condone it but this cleans up the warnings
    warnings.filterwarnings("ignore")
    
    
    #Fix input formats:
    rise_time_depths=[rise_time_depths0,rise_time_depths1]
    #Load the source
    mpi_rupt=home+project_name+'/output/ruptures/mpi_rupt.'+str(rank)+'.'+rupture_name
    fault=genfromtxt(mpi_rupt)  
    
    #Onset times for each subfault
    onset_times=fault[:,12]
    
    #load velocity structure
    structure=genfromtxt(home+project_name+'/structure/'+model_name)
    
    #Frequencies vector
    f=logspace(log10(1/total_duration),log10(1/(2*hf_dt))+0.01,100)
    omega=2*pi*f
    
    #Output time vector (0 is origin time)
    t=arange(0,total_duration,hf_dt)
    
    #Projection object for distance calculations
    g=Geod(ellps='WGS84')
    
    #Create taup velocity model object, paste on top of iaspei91
    #taup_create.build_taup_model(home+project_name+'/structure/bbp_norcal.tvel',output_folder=home+project_name+'/structure/')
#    velmod=TauPyModel(model=home+project_name+'/structure/iquique',verbose=True)
    velmod = TauPyModel(model=home+project_name+'/structure/'+model_name.split('.')[0])
    
    #Get epicentral time
    epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
    
    #Moments
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    subfault_M0=slip*fault[:,10]*fault[:,11]*fault[:,13]
    subfault_M0=subfault_M0*1e7 #to dyne-cm
    relative_subfault_M0=subfault_M0/M0
    Mw=(2./3)*(log10(M0*1e-7)-9.1)
    
    #Corner frequency scaling
    i=where(slip>0)[0] #Non-zero faults
    dl=mean((fault[:,10]+fault[:,11])/2) #predominant length scale
    dl=dl/1000 # to km
    
    #Tau=p perturbation
    tau_perturb=0.1
    
    #Deep faults receive a higher stress
    stress_multiplier=1

    #initalize output seismogram
    tr=Trace()
    tr.stats.station=sta
    tr.stats.delta=hf_dt
    tr.stats.starttime=time_epi
    #info for sac header (added at the end)
    az,backaz,dist_m=g.inv(epicenter[0],epicenter[1],sta_lon,sta_lat)
    dist_in_km=dist_m/1000.    
    
    hf=zeros(len(t))
    
    #Loop over subfaults
    for kfault in range(len(fault)):
        if rank==0:
            #Print status to screen            
            if kfault % 25 == 0:
                if kfault==0:
                    stdout.write('      [.')
                    stdout.flush()
                stdout.write('.')
                stdout.flush()
            if kfault==len(fault)-1:
                stdout.write('.]\n')
                stdout.flush()                
        
        #Include only subfaults with non-zero slip
        if subfault_M0[kfault]>0:
            
            #Get subfault to station distance
            lon_source=fault[kfault,1]
            lat_source=fault[kfault,2]
            azimuth,baz,dist=g.inv(lon_source,lat_source,sta_lon,sta_lat)
            dist_in_degs=kilometer2degrees(dist/1000.)
            
            #Source depth?
            z_source=fault[kfault,3]
            
            #No change
            stress=stress_parameter
            
            #Is subfault in an SMGA?
            #SMGA1
#            radius_in_km=15.0
#            smga_center_lon=-71.501
#            smga_center_lat=-30.918
            
            
            #SMGA2
#            radius_in_km=15.0
#            smga_center_lon=-71.863
#            smga_center_lat=-30.759
            
            #smga3
#            radius_in_km=7.5  
#            smga_center_lon=-72.3923
#            smga_center_lat=-30.58
            
            
            #smga4
            # radius_in_km=7.5  
            # smga_center_lon=-72.3923
            # smga_center_lat=-30.61
            
            # in_smga=is_subfault_in_smga(lon_source,lat_source,smga_center_lon,smga_center_lat,radius_in_km)
            
            # ###Apply multiplier?
            # if in_smga==True:
            #     stress=stress_parameter*stress_multiplier
            #     print("%.4f,%.4f is in SMGA, stress is %d" % (lon_source,lat_source,stress))
            # else:
            #     stress=stress_parameter
            
            #Apply multiplier?
            #if slip[kfault]>7.5:
            #    stress=stress_parameter*stress_multiplier
            ##elif lon_source>-72.057 and lon_source<-71.2 and lat_source>-30.28:
            ##    stress=stress_parameter*stress_multiplier
            #else:
            #    stress=stress_parameter
                
            #Apply multiplier?
            #if z_source>high_stress_depth:
            #    stress=stress_parameter*stress_multiplier
            #else:
            #    stress=stress_parameter
            
            # Frankel 95 scaling of corner frequency #verified this looks the same in GP
            # Right now this applies the same factor to all faults
            fc_scale=(M0)/(N*stress*dl**3*1e21) #Frankel scaling
            small_event_M0 = stress*dl**3*1e21
            
        

            
            #Get rho, alpha, beta at subfault depth
            zs=fault[kfault,3]
            mu,alpha,beta=get_mu(structure,zs,return_speeds=True)
            rho=mu/beta**2
            
            #Get radiation scale factor
            Spartition=1/2**0.5
            if component=='N' :
                component_angle=0
            elif component=='E':
                component_angle=90
            
            rho=rho/1000 #to g/cm**3
            beta=(beta/1000)*1e5 #to cm/s
            alpha=(alpha/1000)*1e5
            
            # print('rho = '+str(rho))
            # print('beta = '+str(beta))
            # print('alpha = '+str(alpha))
            
            #Verified this produces same value as in GP
            CS=(2*Spartition)/(4*pi*(rho)*(beta**3))
            CP=2/(4*pi*(rho)*(alpha**3))

            
            #Get local subfault rupture speed
            beta=beta/100 #to m/s
            vr=hfsims.get_local_rupture_speed(zs,beta,rise_time_depths)
            vr=vr/1000 #to km/s
            dip_factor=hfsims.get_dip_factor(fault[kfault,5],fault[kfault,8],fault[kfault,9])
            
            #Subfault corner frequency
            c0=2.0 #GP2015 value
            fc_subfault=(c0*vr)/(dip_factor*pi*dl)
            
            #get subfault source spectrum
            #S=((relative_subfault_M0[kfault]*M0/N)*f**2)/(1+fc_scale*(f/fc_subfault)**2)
            S=small_event_M0*(omega**2/(1+(f/fc_subfault)**2))
            frankel_conv_operator= fc_scale*((fc_subfault**2+f**2)/(fc_subfault**2+fc_scale*f**2))
            S=S*frankel_conv_operator
            
            #get high frequency decay
            P=exp(-pi*kappa*f)
            
            
            #Get other geometric parameters necessar for radiation pattern
            strike=fault[kfault,4]
            dip=fault[kfault,5]
            ss=fault[kfault,8]
            ds=fault[kfault,9]
            rake=rad2deg(arctan2(ds,ss))
            
            #Get ray paths for all direct P arrivals
            Ppaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['P','p'])
            
            #Get ray paths for all direct S arrivals
            try:
                Spaths=velmod.get_ray_paths(zs,dist_in_degs,phase_list=['S','s'])
            except:
                Spaths=velmod.get_ray_paths(zs+tau_perturb,dist_in_degs,phase_list=['S','s'])
            
            #sometimes there's no S, weird I know. Check twice.
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+tau_perturb,dist_in_degs,phase_list=['S','s'])
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+5*tau_perturb,dist_in_degs,phase_list=['S','s'])   
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs-5*tau_perturb,dist_in_degs,phase_list=['S','s'])
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+5*tau_perturb,dist_in_degs,phase_list=['S','s'])  
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs-10*tau_perturb,dist_in_degs,phase_list=['S','s'])
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+10*tau_perturb,dist_in_degs,phase_list=['S','s']) 
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs-50*tau_perturb,dist_in_degs,phase_list=['S','s'])
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+50*tau_perturb,dist_in_degs,phase_list=['S','s']) 
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs-75*tau_perturb,dist_in_degs,phase_list=['S','s'])
            if len(Spaths)==0:
                Spaths=velmod.get_ray_paths(zs+75*tau_perturb,dist_in_degs,phase_list=['S','s']) 
            if len(Spaths)==0:
                print('ERROR: I give up, no direct S in spite of multiple attempts at subfault '+str(kfault))

            #Get direct s path and moho reflection
            mohoS=None
            directS=Spaths[0]
            directP=Ppaths[0]
            #print len(Spaths)
            if len(Spaths)==1: #only direct S
                pass
            else:
                #turn_depth=zeros(len(Spaths)-1) #turning depth of other non-direct rays
                #for k in range(1,len(Spaths)):
                #    turn_depth[k-1]=Spaths[k].path['depth'].max()
                ##If there's a ray that turns within 2km of Moho, callt hat guy the Moho reflection
                #deltaz=abs(turn_depth-moho_depth_in_km)
                #i=argmin(deltaz)
                #if deltaz[i]<2: #Yes, this is a moho reflection
                #    mohoS=Spaths[i+1]
                #else:
                #    mohoS=None
                mohoS=Spaths[-1]
                

                 
 
            #######         Build Direct P ray           ######
            if Pwave==True:
                take_off_angle_P=directP.takeoff_angle
                
                # #Get attenuation due to geometrical spreading (from the path length)
                # path_length_P=hfsims.get_path_length(directP,zs,dist_in_degs)
                # path_length_P=path_length_P*100 #to cm
                
                # #Get effect of intrinsic attenuation for that ray (path integrated)
                # #Q_P=hfsims.get_attenuation(f,structure,directP,Qexp,Qtype='P')   <- This causes problems and I don't know why underlying assumptions might be bad
                # Q_P=hfsims.get_attenuation(f,structure,directS,Qexp,Qtype='S')
                
                # #get quarter wavelength amplificationf actors
                # # pass rho in kg/m^3 (this units nightmare is what I get for following Graves' code)
                # I_P=hfsims.get_amplification_factors(f,structure,zs,alpha,rho*1000)
                
                # #Build the entire path term
                # G_P=(I_P*Q_P)/path_length_P
                
                #Get attenuation due to geometrical spreading (from the path length)
                path_length_S=hfsims.get_path_length(directS,zs,dist_in_degs)
                path_length_S=path_length_S*100 #to cm
                
                #Get effect of intrinsic aptimeenuation for that ray (path integrated)
                Q_S=hfsims.get_attenuation(f,structure,directS,Qexp)
                
                #get quarter wavelength amplificationf actors
                # pass rho in kg/m^3 (this units nightmare is what I get for following Graves' code)
                I_S=hfsims.get_amplification_factors(f,structure,zs,beta,rho*1000)
                
                #Build the entire path term
                G_S=(I_S*Q_S)/path_length_S
                
                

                #Get conically averaged radiation pattern terms
                RP=hfsims.conically_avg_P_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_P)
                RP=abs(RP)
                   
                #Get partition of Pwave into Z and N,E components 
                incidence_angle=directP.incident_angle
                Npartition,Epartition,Zpartition=hfsims.get_P_wave_partition(incidence_angle,azimuth)
                if component=='Z':
                   Ppartition=Zpartition 
                elif component=='N':
                    Ppartition=Npartition
                else:
                    Ppartition=Epartition
                    
                #And finally multiply everything together to get the subfault amplitude spectrum
                AP=CP*S*G_S*P*RP*Ppartition           

                #Generate windowed time series
                duration=1./fc_subfault+0.09*(dist/1000)
                w=hfsims.windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
                
                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
                hf_seis_P=hfsims.apply_spectrum(w,AP,f,hf_dt)
                
                #save thigns to check
                # if sta=='AL2H':
                #     path_out = '/Users/dmelgarm/FakeQuakes/ONC_debug/analysis/frequency/Pwave/'
                #     path_out = path_out+str(kfault)
                #     # savetxt(path_out+'.all',c_[f,AP])
                #     # savetxt(path_out+'.source',c_[f,CP*S])
                #     # savetxt(path_out+'.path',c_[f,G_P])
                #     # savetxt(path_out+'.site',c_[f,P])
                
                
                #What time after OT should this time series start at?
                time_insert=directP.path['time'][-1]+onset_times[kfault]
                i=argmin(abs(t-time_insert))
                j=i+len(hf_seis_P)
                
                #Check seismogram doesn't go past last sample
                if i<len(hf)-1: #if i (the beginning of the seimogram) is less than the length
                    if j>len(hf): #seismogram goes past total_duration length, trim it
                        len_paste=len(hf)-i
                        j=len(hf)
                        #Add seismogram
                        hf[i:j]=hf[i:j]+real(hf_seis_P[0:len_paste])
                    else: #Lengths are fine
                        hf[i:j]=hf[i:j]+real(hf_seis_P)      
                else: #Seismogram starts after end of available space
                    pass   
                
                                           
                                                                  
                                                                                                                
                          
            #######         Build Direct S ray           ######

            if Swave==True:
                take_off_angle_S=directS.takeoff_angle
                
                #Get attenuation due to geometrical spreading (from the path length)
                path_length_S=hfsims.get_path_length(directS,zs,dist_in_degs)
                path_length_S=path_length_S*100 #to cm
                
                #Get effect of intrinsic aptimeenuation for that ray (path integrated)
                Q_S=hfsims.get_attenuation(f,structure,directS,Qexp)
                
                #get quarter wavelength amplificationf actors
                # pass rho in kg/m^3 (this units nightmare is what I get for following Graves' code)
                I_S=hfsims.get_amplification_factors(f,structure,zs,beta,rho*1000)
                
                #Build the entire path term
                G_S=(I_S*Q_S)/path_length_S
    
                #Get conically averaged radiation pattern terms
                if component=='Z':
                    RP_vert=hfsims.conically_avg_vert_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_S)
                    #And finally multiply everything together to get the subfault amplitude spectrum
                    AS=CS*S*G_S*P*RP_vert   
                    # print('... RP_vert = '+str(RP_vert))
                else:
                    RP=hfsims.conically_avg_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_S,component_angle)
                    RP=abs(RP)
                    # print('... RP_horiz = '+str(RP))
                    #And finally multiply everything together to get the subfault amplitude spectrum
                    AS=CS*S*G_S*P*RP                
    
                #Generate windowed time series
                duration=1./fc_subfault+0.063*(dist/1000)
                w=hfsims.windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
                #w=windowed_gaussian(3*duration,hf_dt,window_type='cua',ptime=Ppaths[0].path['time'][-1],stime=Spaths[0].path['time'][-1])
                
                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
                hf_seis_S=hfsims.apply_spectrum(w,AS,f,hf_dt)
                
                #save thigns to check
                # if sta=='AL2H':
                #     path_out = '/Users/dmelgarm/FakeQuakes/ONC_debug/analysis/frequency/Swave/'
                #     path_out = path_out+str(kfault)
                #     # savetxt(path_out+'.soverp',c_[f,(CS*S)/(CP*S)])
                #     savetxt(path_out+'.all',c_[f,AS])
                #     savetxt(path_out+'.source',c_[f,CS*S])
                #     savetxt(path_out+'.path',c_[f,G_S])
                #     savetxt(path_out+'.site',c_[f,P])
                
                #What time after OT should this time series start at?
                time_insert=directS.path['time'][-1]+onset_times[kfault]
                #print 'ts = '+str(time_insert)+' , Td = '+str(duration)
                #time_insert=Ppaths[0].path['time'][-1]
                i=argmin(abs(t-time_insert))
                j=i+len(hf_seis_S)
                
                
                #Check seismogram doesn't go past last sample
                if i<len(hf)-1: #if i (the beginning of the seimogram) is less than the length
                    if j>len(hf): #seismogram goes past total_duration length, trim it
                        len_paste=len(hf)-i
                        j=len(hf)
                        #Add seismogram
                        hf[i:j]=hf[i:j]+real(hf_seis_S[0:len_paste])
                    else: #Lengths are fine
                        hf[i:j]=hf[i:j]+real(hf_seis_S)
                else: #Beginning of seismogram is past end of available space
                    pass

            
            
            #######         Build Moho reflected S ray           ######
#            if mohoS==None:
#                pass
#            else:
#                if kfault%100==0:
#                    print '... ... building Moho reflected S wave'
#                take_off_angle_mS=mohoS.takeoff_angle
#                
#                #Get attenuation due to geometrical spreading (from the path length)
#                path_length_mS=get_path_length(mohoS,zs,dist_in_degs)
#                path_length_mS=path_length_mS*100 #to cm
#                
#                #Get effect of intrinsic aptimeenuation for that ray (path integrated)
#                Q_mS=get_attenuation(f,structure,mohoS,Qexp)
#                
#                #Build the entire path term
#                G_mS=(I*Q_mS)/path_length_mS
#
#                #Get conically averaged radiation pattern terms
#                if component=='Z':
#                    RP_vert=conically_avg_vert_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_mS)
#                    #And finally multiply everything together to get the subfault amplitude spectrum
#                    A=C*S*G_mS*P*RP_vert   
#                else:
#                    RP=conically_avg_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_mS,component_angle)
#                    RP=abs(RP)
#                    #And finally multiply everything together to get the subfault amplitude spectrum
#                    A=C*S*G_mS*P*RP                
#
#                #Generate windowed time series
#                duration=1./fc_subfault+0.063*(dist/1000)
#                w=windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
#                #w=windowed_gaussian(3*duration,hf_dt,window_type='cua',ptime=Ppaths[0].path['time'][-1],stime=Spaths[0].path['time'][-1])
#                
#                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
#                hf_seis=apply_spectrum(w,A,f,hf_dt)
#                
#                #What time after OT should this time series start at?
#                time_insert=mohoS.path['time'][-1]+onset_times[kfault]
#                #print 'ts = '+str(time_insert)+' , Td = '+str(duration)
#                #time_insert=Ppaths[0].path['time'][-1]
#                i=argmin(abs(t-time_insert))
#                j=i+len(hf_seis)
#                
#                #Add seismogram
#                hf[i:j]=hf[i:j]+hf_seis
#                
#                #Done, reset
#                mohoS=None        
           
    #Done
    tr.data=hf/100 #convert to m/s**2
    #Add station location, event location, and first P-wave arrival time to SAC header
    tr.stats.update({'sac':{'stlo':sta_lon,'stla':sta_lat,'evlo':epicenter[0],'evla':epicenter[1],'evdp':epicenter[2],'dist':dist_in_km,'az':az,'baz':backaz,'mag':Mw}}) #,'idep':"ACC (m/s^2)" not sure why idep won't work
    
    #Write out to file 
    #old
    rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new
    rupture=rupture_name.rsplit('.',1)[0]
    if not path.exists(home+project_name+'/output/waveforms/'+rupture+'/'):
        makedirs(home+project_name+'/output/waveforms/'+rupture+'/')
    if rank < 10:
        tr.write(home+project_name+'/output/waveforms/'+rupture+'/'+sta+'.HN'+component+'.00'+str(rank)+'.sac',format='SAC')
    elif rank < 100:
        tr.write(home+project_name+'/output/waveforms/'+rupture+'/'+sta+'.HN'+component+'.0'+str(rank)+'.sac',format='SAC')
    else:
        tr.write(home+project_name+'/output/waveforms/'+rupture+'/'+sta+'.HN'+component+'.'+str(rank)+'.sac',format='SAC')

    
#If main entry point
if __name__ == '__main__':
    import sys
    from mpi4py import MPI
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #Map command line arguments to funciton arguments:
    if sys.argv[1]=='run_parallel_hfsims':
        #Parse command line arguments
        home=sys.argv[2]
        project_name=sys.argv[3]
        rupture_name=sys.argv[4]
        N=int(sys.argv[5])
        M0=float(sys.argv[6])
        sta=sys.argv[7]
        sta_lon=float(sys.argv[8])
        sta_lat=float(sys.argv[9])
        model_name=sys.argv[10]
        rise_time_depths0=int(sys.argv[11])
        rise_time_depths1=int(sys.argv[12])
        moho_depth_in_km=float(sys.argv[13])
        component=sys.argv[14]
        total_duration=int(sys.argv[15])
        hf_dt=float(sys.argv[16])
        stress_parameter=float(sys.argv[17])
        kappa=float(sys.argv[18])
        Qexp=float(sys.argv[19])
        Pwave=sys.argv[20]
        if Pwave=='True':
            Pwave=True
        Swave=sys.argv[21]
        if Swave=='True':
            Swave=True
        high_stress_depth=int(sys.argv[22])
        run_parallel_hfsims(home,project_name,rupture_name,N,M0,sta,sta_lon,sta_lat,component,model_name,rise_time_depths0,rise_time_depths1,moho_depth_in_km,total_duration,hf_dt,stress_parameter,kappa,Qexp,Pwave,Swave,high_stress_depth,rank,size)
    else:
        print("ERROR: You're not allowed to run "+sys.argv[1]+" from the shell or it does not exist")
        