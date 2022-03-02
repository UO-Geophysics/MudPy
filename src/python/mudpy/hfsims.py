def stochastic_simulation(home,project_name,rupture_name,sta,sta_lon,sta_lat,component,model_name,
        rise_time_depths,moho_depth_in_km,total_duration=100,hf_dt=0.01,stress_parameter=50,
        kappa=0.04,Qexp=0.6,Pwave=False,Swave=True,high_stress_depth=1e4): 
    '''
    Run stochastic HF sims
    
    stress parameter is in bars
    '''
    
    from numpy import genfromtxt,pi,logspace,log10,mean,where,exp,arange,zeros,argmin,rad2deg,arctan2,real
    from pyproj import Geod
    from obspy.geodetics import kilometer2degrees
    from obspy.taup import TauPyModel
    from mudpy.forward import get_mu, write_fakequakes_hf_waveforms_one_by_one,read_fakequakes_hypo_time
    from obspy import Stream,Trace
    from sys import stdout
    import warnings


    #print out what's going on:
    out='''Running with input parameters:
    home = %s
    project_name = %s
    rupture_name = %s
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
    '''%(home,project_name,rupture_name,sta,str(sta_lon),str(sta_lat),model_name,str(rise_time_depths),
    str(moho_depth_in_km),str(total_duration),str(hf_dt),str(stress_parameter),
    str(kappa),str(Qexp),str(component),str(Pwave),str(Swave),str(high_stress_depth))
    print(out)

#    rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
#    log=home+project_name+'/output/waveforms/'+rupture+'/'+sta+'.HN'+component+'.1cpu.log'
#    logfile=open(log,'w')
#    logfile.write(out)
    #print 'stress is '+str(stress_parameter)

    #I don't condone it but this cleans up the warnings
    warnings.filterwarnings("ignore")
    
    #Load the source
    fault=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)    
    
    #Onset times for each subfault
    onset_times=fault[:,12]
    
    #load velocity structure
    structure=genfromtxt(home+project_name+'/structure/'+model_name)
    
    #Frequencies vector
    f=logspace(log10(hf_dt),log10(1/(2*hf_dt))+0.01,50)
    omega=2*pi*f
    
    #Output time vector (0 is origin time)
    t=arange(0,total_duration,hf_dt)
    
    #Projection object for distance calculations
    g=Geod(ellps='WGS84')
    
    #Create taup velocity model object, paste on top of iaspei91
    #taup_create.build_taup_model(home+project_name+'/structure/bbp_norcal.tvel',output_folder=home+project_name+'/structure/')
    velmod=TauPyModel(model=home+project_name+'/structure/iquique',verbose=True)
    #Get epicentral time
    epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
    
    #Moments
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    subfault_M0=slip*fault[:,10]*fault[:,11]*fault[:,13]
    subfault_M0=subfault_M0*1e7 #to dyne-cm
    M0=subfault_M0.sum()
    relative_subfault_M0=subfault_M0/M0
    Mw=(2./3)*(log10(M0*1e-7)-9.1)
    
    #Corner frequency scaling
    i=where(slip>0)[0] #Non-zero faults
    N=len(i) #number of subfaults
    dl=mean((fault[:,10]+fault[:,11])/2) #predominant length scale
    dl=dl/1000 # to km
    
    #Tau=p perturbation
    tau_perturb=0.1
    
    #Deep faults receive a higher stress
    stress_multiplier=3
         
    print('... working on '+component+' component semistochastic waveform for station '+sta)

    #initalize output seismogram
    tr=Trace()
    tr.stats.station=sta
    tr.stats.delta=hf_dt
    tr.stats.starttime=time_epi
    #info for sac header (added at the end)
    az,backaz,dist_m=g.inv(epicenter[0],epicenter[1],sta_lon,sta_lat)
    dist_in_km=dist_m/1000.    
    
    hf=zeros(len(t))
    
#    out='''Parameters before we get into subfault calculations:
#    rupture_name = %s
#    epicenter = %s
#    time_epi = %s
#    M0 = %E
#    Mw = %10.4f
#    Num_Subfaults = %i
#    dl = %.2f
#    Dist_in_km = %10.4f
#    '''%(rupture_name,str(epicenter),str(time_epi),M0,Mw,int(N),dl,dist_in_km)
#    print out
#    logfile.write(out)
    
    #Loop over subfaults
#    earliestP=1e10  #something outrageously high
#    earliestP_kfault=1e10
    for kfault in range(len(fault)):
        
        #Print status to screen            
        if kfault % 150 == 0:
            if kfault==0:
                stdout.write('      [')
                stdout.flush()
            stdout.write('.')
            stdout.flush()
        if kfault==len(fault)-1:
            stdout.write(']\n')
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
            #radius_in_km=15.0
            #smga_center_lon=-69.709200
            #smga_center_lat=-19.683600
            #in_smga=is_subfault_in_smga(lon_source,lat_source,smga_center_lon,smga_center_lat,radius_in_km)
            #
            ###Apply multiplier?
            #if in_smga==True:
            #    stress=stress_parameter*stress_multiplier
            #    print "%.4f,%.4f is in SMGA, stress is %d" % (lon_source,lat_source,stress)
            #else:
            #    stress=stress_parameter
            
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
            
            #Verified this produces same value as in GP
            CS=(2*Spartition)/(4*pi*(rho)*(beta**3))
            CP=2/(4*pi*(rho)*(alpha**3))
            
            #Get local subfault rupture speed
            beta=beta/100 #to m/s
            vr=get_local_rupture_speed(zs,beta,rise_time_depths)
            vr=vr/1000 #to km/s
            dip_factor=get_dip_factor(fault[kfault,5],fault[kfault,8],fault[kfault,9])
            
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
            
            
#            if kfault==0:
#                out='''Parameters within subfault calculations:
#                kfault_lon = %10.4f
#                kfault_lat = %10.4f
#                CS = %s
#                CP = %s
#                S[0] = %s
#                frankel_conv_operator[0] = %s
#                '''%(fault[kfault,1],fault[kfault,2],str(CS),str(CP),str(S[0]),str(frankel_conv_operator[0]))
#                print out
#                logfile.write(out)
            
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
                
                #Get attenuation due to geometrical spreading (from the path length)
                path_length_P=get_path_length(directP,zs,dist_in_degs)
                path_length_P=path_length_P*100 #to cm
                
                #Get effect of intrinsic attenuation for that ray (path integrated)
                Q_P=get_attenuation(f,structure,directP,Qexp,Qtype='P')
                
                #get quarter wavelength amplificationf actors
                # pass rho in kg/m^3 (this units nightmare is what I get for following Graves' code)
                I_P=get_amplification_factors(f,structure,zs,alpha,rho*1000)
            
                #Build the entire path term
                G_P=(I_P*Q_P)/path_length_P

                #Get conically averaged radiation pattern terms
                RP=conically_avg_P_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_P)
                RP=abs(RP)
                   
                #Get partition of Pwave into Z and N,E components 
                incidence_angle=directP.incident_angle
                Npartition,Epartition,Zpartition=get_P_wave_partition(incidence_angle,azimuth)
                if component=='Z':
                   Ppartition=Zpartition 
                elif component=='N':
                    Ppartition=Npartition
                else:
                    Ppartition=Epartition
                    
                #And finally multiply everything together to get the subfault amplitude spectrum
                AP=CP*S*G_P*P*RP*Ppartition           

                #Generate windowed time series
                duration=1./fc_subfault+0.09*(dist/1000)
                w=windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
                
                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
                hf_seis_P=apply_spectrum(w,AP,f,hf_dt)
                
                #What time after OT should this time series start at?
                time_insert=directP.path['time'][-1]+onset_times[kfault]
#                if directP.time+onset_times[kfault] < earliestP:
#                    earliestP=directP.time+onset_times[kfault]
#                    earliestP_kfault=kfault
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
                path_length_S=get_path_length(directS,zs,dist_in_degs)
                path_length_S=path_length_S*100 #to cm
                
                #Get effect of intrinsic aptimeenuation for that ray (path integrated)
                Q_S=get_attenuation(f,structure,directS,Qexp)
                
                #get quarter wavelength amplificationf actors
                # pass rho in kg/m^3 (this units nightmare is what I get for following Graves' code)
                I_S=get_amplification_factors(f,structure,zs,beta,rho*1000)
                #Build the entire path term
                G_S=(I_S*Q_S)/path_length_S
    
                #Get conically averaged radiation pattern terms
                if component=='Z':
                    RP_vert=conically_avg_vert_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_S)
                    print(RP_vert)
                    #And finally multiply everything together to get the subfault amplitude spectrum
                    AS=CS*S*G_S*P*RP_vert   
                else:
                    RP=conically_avg_radiation_pattern(strike,dip,rake,azimuth,take_off_angle_S,component_angle)
                    RP=abs(RP)
                    print(RP)
                    #And finally multiply everything together to get the subfault amplitude spectrum
                    AS=CS*S*G_S*P*RP                
    
                #Generate windowed time series
                duration=1./fc_subfault+0.063*(dist/1000)
                w=windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
                #w=windowed_gaussian(3*duration,hf_dt,window_type='cua',ptime=Ppaths[0].path['time'][-1],stime=Spaths[0].path['time'][-1])
                
                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
                hf_seis_S=apply_spectrum(w,AS,f,hf_dt)
                
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
#        if kfault==0:
#            out=''' More:
#            fc_scale = %10.4f
#            subfaultM0 = %E
#            mu = %E
#            CS = %E
#            CP = %E
#            vr = %10.4f
#            dip_factor = %10.4f
#            fc_subfault = %10.4f
#            directS = %s
#            directP = %s
#            '''%(fc_scale,subfault_M0[kfault],mu,CS,CP,vr,dip_factor,fc_subfault,str(directS.time),str(directP.time))
#            print out
#            logfile.write(out)
#    logfile.close()
    #Done
    tr.data=hf/100 #convert to m/s**2
    #Add station location, event location, and first P-wave arrival time to SAC header
    tr.stats.update({'sac':{'stlo':sta_lon,'stla':sta_lat,'evlo':epicenter[0],'evla':epicenter[1],'evdp':epicenter[2],'dist':dist_in_km,'az':az,'baz':backaz,'mag':Mw}}) #,'idep':"ACC (m/s^2)" not sure why idep won't work
    #Return trace for writing to file    
#    print "Earliest P wave Comes at " + str(earliestP) + "after OT, from location " + str(fault[earliestP_kfault,1]) + ", " + str(fault[earliestP_kfault,2]) + ", " +str(fault[earliestP_kfault,3])
    return tr
            
                        
                        
def get_local_rupture_speed(zs,beta,rise_time_depths): 
    '''
    Get local rupture speed
    '''
    
    if zs<rise_time_depths[0]:
        vr=0.56*beta
    elif zs>rise_time_depths[1]:
        vr=0.8*beta
    else:
        m=(0.24/3.0)
        b=0.8-8*m
        multiplier=m*zs+b
        vr=multiplier*beta
    return vr


def get_dip_factor(dip,ss,ds):
    '''
    This now uses GP2015
    '''
    
    from numpy import arctan2,rad2deg
    
    if dip>45:
        FD=1-(dip-45)/45.
    else:
        FD=1
    
    rake=rad2deg(arctan2(ds,ss))
    
    if (rake>0 and rake <180):
        FR=1-(rake-90)/90.
    else:
        FR=0
        
    dip_factor=1./(1+FD*FR*0.1)

    return dip_factor                    
                                                                                                                        

                                                            
def get_amplification_factors(f,structure,zs,beta,rho):
    '''
    Get quarter wavelength amplification factors, this guy operates in SI units
    '''

    from numpy import zeros,arange
    from scipy.interpolate import interp1d

    #get mean velocity to all depths
    zmax=200.0*1000 #in m
    z=zeros(len(structure)*2)
    vs=zeros(len(structure)*2)
    rho_model=zeros(len(structure)*2)
    for k in range(1,len(structure)):
        z[2*k-1]=z[2*k-2]+structure[k-1,0]
        z[2*k]=z[2*k-1]+0.000001
        vs[2*k-2]=structure[k-1,1]
        vs[2*k-1]=structure[k-1,1]
        rho_model[2*k-2]=structure[k-1,3]
        rho_model[2*k-1]=structure[k-1,3]
        
    z[-2]=z[-3]+0.00001
    z=z*1000
    z[-1]=zmax
    
    vs[-1]=structure[-1,1]
    vs[-2]=structure[-1,1]
    rho_model[-1]=structure[-1,3]
    rho_model[-2]=structure[-1,3]
    
    #interpolate
    interpolator=interp1d(z,vs,bounds_error=False)
    interpolator2=interp1d(z,rho_model,bounds_error=False)
    dz=1.0
    zinterp=arange(0,zmax,dz)
    vsinterp=interpolator(zinterp)
    rhointerp=interpolator2(zinterp)
    
    #mean velocity
    mean_vs=vsinterp.cumsum()/arange(1,len(vsinterp)+1)*1000
    
    #mean rho
    mean_rho=rhointerp.cumsum()/arange(1,len(rhointerp)+1)*1000
    
    #frequency for each depth
    fz=mean_vs/(4*zinterp)

    #amplifications at those frequencies
    Afz=((beta*rho)/(mean_rho*mean_vs))**0.5
    
    #resample to frequencies of interest
    interpolator=interp1d(fz,Afz,bounds_error=False)
    I=interpolator(f)
    
    return I
        

def get_path_length(ray,zs,dist_in_degs):
    
    from numpy import diff
    
    radius_of_earth=6371e3
    dist=ray.path['dist']
    dist=dist*radius_of_earth    #now this is in meters
    depth=ray.path['depth']*1000 #this is to in meters now
    path_dist=(diff(dist)**2+diff(depth)**2)**0.5
    path_length=path_dist.sum()
    
    return path_length
    

def get_attenuation(f,structure,ray,Qexp,Qtype='S'):
    '''
    Get effect of intrinsic aptimeenuation along the ray path
    '''
    
    from numpy import diff,zeros,exp,pi
    from mudpy.forward import get_Q
    
    time=ray.path['time']
    time_in_layer=diff(time)
    depth=ray.path['depth']
    omega=2*pi*f
    
    Qp=zeros(len(time_in_layer))
    Qs=zeros(len(time_in_layer))
    
    for k in range(len(Qp)):
        Qp[k],Qs[k]=get_Q(structure,depth[k])

        
    #Get the travel tiem weighted sum
    if Qtype=='S':
        weightedQ=sum(time_in_layer/Qs)
        
    else:
        weightedQ=sum(time_in_layer/Qp)

    
    #get frequency dependence
    Q=exp(-pi*weightedQ*f**(1-Qexp))
    
    return Q


def windowed_gaussian(duration,hf_dt,window_type='saragoni_hart',M=5.0,dist_in_km=50.,std=1.0,ptime=10,stime=20):
    '''
    Get a gaussian white noise time series and window it
    '''
    
    from numpy.random import normal
    from numpy import log,exp,arange
    from scipy.special import gamma
    
    mean=0.0
    num_samples = int(duration/hf_dt)
    #If num_smaples is even then make odd for FFT stuff later on
    if num_samples%2==0:
        num_samples+=1
    t=arange(0,num_samples*hf_dt,hf_dt)
    t=t[0:num_samples]

    noise = normal(mean, std, size=num_samples)
    
    if window_type=='saragoni_hart':
        epsilon=0.2
        eta=0.05
        b=-epsilon*log(eta)/(1+eta*(log(epsilon)-1))
        c=b/(epsilon*duration)
        a=(((2*c)**(2*b+1))/gamma(2*b+1))**0.5
        window=a*t**b*exp(-c*t)
    elif window_type=='cua':
        ptime=0
        window=cua_envelope(M,dist_in_km,t,ptime,stime,Pcoeff=0,Scoeff=12)
    elif window_type==None: #jsut white noise, no window
        window=1
        
        
    noise=noise*window
    
    return noise

        
def apply_spectrum(w,A,f,hf_dt,is_gnss=False,gnss_scale=1/2**0.5):
    '''
    Apply the modeled spectrum to the windowed time series
    
    '''
    
    from numpy import fft,angle,cos,sin,mean,zeros,real
    from scipy.interpolate import interp1d
    
    #to frequency domain
    fourier=fft.fft(w)
    freq=fft.fftfreq(len(w),hf_dt)
    
    #Get positive frequencies
    Nf=len(freq)
    positive_freq=freq[1:int(1+Nf/2)]
      
    #Make POWER spectrum of windowed time series have a mean of 1
    #norm_factor=hf_dt*mean(abs(fourier)**2)**0.5
    #norm_factor=mean(abs(fourier))
    norm_factor=mean(abs(fourier)**2)**0.5
    fourier=fourier/norm_factor
    
    #Keep phase
    phase=angle(fourier)
    
    #resample model amplitude spectr to frequencies
    interp=interp1d(f,A,bounds_error=False)
    amplitude_positive=interp(positive_freq)
    
    #Place in correct order A[0] is DC value then icnreasing positive freq then decreasing negative freq
    amplitude=zeros(len(freq))
    #DC value
    amplitude[0]=0
    #Positive freq. div by 2 to keep right power
    amplitude[1:int(1+Nf/2)]=amplitude_positive/2
    #Negative freq
    amplitude[int(1+Nf/2):]=amplitude_positive[::-1]/2
    
    #Scale for GNSS displacememnts?
    if is_gnss:
        amplitude /= gnss_scale 
    
    #Apply model amplitude spectrum
    amplitude=amplitude*abs(fourier)
    
    #Obtain complex foureier series
    R=amplitude*cos(phase)
    I=amplitude*sin(phase)
    fourier=R+I*1j
    
    #ifft
    seis=real(fft.ifft(fourier))
    
    if is_gnss:
        seis *= len(seis)**0.5
    else:
        seis=seis*len(seis)
    
    return seis         
                                
    
 
def cua_envelope(M,dist_in_km,times,ptime,stime,Pcoeff=0,Scoeff=12):
    '''
    Cua envelopes, modified from Ran Nof's Cua2008 module
    '''
    from numpy import where,sqrt,exp,log10,arctan,pi,zeros
    
    a = [0.719, 0.737, 0.801, 0.836, 0.950, 0.943, 0.745, 0.739, 0.821, 0.812, 0.956, 0.933,
            0.779, 0.836, 0.894, 0.960, 1.031, 1.081, 0.778, 0.751, 0.900, 0.882, 1.042, 1.034]
    b = [-3.273e-3, -2.520e-3, -8.397e-4, -5.409e-4, -1.685e-6, -5.171e-7, -4.010e-3, -4.134e-3,
                -8.543e-4, -2.652e-6, -1.975e-6, -1.090e-7, -2.555e-3, -2.324e-3, -4.286e-4, -8.328e-4,
                -1.015e-7, -1.204e-6, -2.66e-5, -2.473e-3, -1.027e-5,- 5.41e-4, -1.124e-5, -4.924e-6]
    d = [-1.195, -1.26, -1.249, -1.284, -1.275, -1.161, -1.200, -1.199, -1.362, -1.483, -1.345, -1.234,
                -1.352, -1.562, -1.440, -1.589, -1.438, -1.556, -1.385, -1.474, -1.505, -1.484, -1.367, -1.363]
    c1 = [1.600, 2.410, 0.761, 1.214, 2.162, 2.266, 1.752, 2.030, 1.148, 1.402, 1.656, 1.515,
                1.478, 2.423, 1.114, 1.982, 1.098, 1.946, 1.763, 1.593, 1.388, 1.530, 1.379, 1.549]
    c2 = [1.045, 0.955, 1.340, 0.978, 1.088, 1.016, 1.091, 1.972, 1.100, 0.995, 1.164, 1.041,
                1.105, 1.054, 1.110, 1.067, 1.133, 1.091, 1.112, 1.106, 1.096, 1.04, 1.178, 1.082]
    e = [-1.065, -1.051, -3.103, -3.135, -4.958, -5.008, -0.955, -0.775, -2.901, -2.551, -4.799, -4.749,
                -0.645, -0.338, -2.602, -2.351, -4.342, -4.101, -0.751, -0.355, -2.778, -2.537, -4.738, -4.569]
    sig_uncorr = [0.307, 0.286, 0.268, 0.263, 0.284, 0.301, 0.288, 0.317, 0.263, 0.298, 02.83, 0.312,
                0.308, 0.312, 0.279, 0.296, 0.277, 0.326, 0.300, 0.300, 0.250, 0.270, 0.253, 0.286]
    sig_corr = [0.233, 0.229, 0.211, 0.219, 0.239, 0.247, 0.243, 0.256, 0.231, 0.239, 0.254, 0.248,
                0.243, 0.248, 0.230, 0.230, 0.233, 0.236, 0.238, 0.235, 0.220, 0.221, 0.232, 0.230]
    
    # Coefficienstime for eqn: log(env_param) = alpha*M + beta*R + delta*logR + mu
    # Coefficienstime and equation for t_rise (rise time):
    
    alpha_t_rise = [0.06, 0.07, 0.06, 0.07, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06, 0.08, 0.067,
                0.064, 0.055, 0.093, 0.087, 0.109, 0.12, 0.069, 0.059, 0.116, 0.11, 0.123, 0.124]  
    beta_t_rise = [5.5e-4, 1.2e-3, 1.33e-3, 4.35e-4, 1.29e-3, 1.19e-3, 7.45e-4, 5.87e-4, 7.32e-4, 1.08e-3, 1.64e-3, 1.21e-3,
                0, 1.21e-3, 0, 4.0e-4, 7.68e-4, 0, 0, 2.18e-3, 0, 1.24e-3, 1.3e-3, 0]
    delta_t_rise = [0.27, 0.24, 0.23, 0.47, 0.27, 0.47, 0.37, 0.23, 0.25, 0.22, 0.13, 0.28,
                0.48, 0.34, 0.48, 0.49, 0.38, 0.45, 0.49, 0.26, 0.503, 0.38, 0.257, 0.439]
    mu_t_rise = [-0.37, -0.38, -0.34, -0.68, -0.34, -0.58, -0.51, -0.37, -0.37, -0.36, -0.33, -0.46,
                -0.89, -0.66, -0.96, -0.98, -0.87,-0.89,-0.97, -0.66, -1.14, -0.91, -0.749, -0.82]
    
    # Coefficienstime and equation for delta_t (wave duration):
    
    alpha_delta_t = [0, 0.03, 0.054, 0.03, 0.047, 0.051, 0, 0, 0.046, 0.031, 0.058, 0.043,
                0, 0.028, 0.02, 0.028, 0.04, 0.03, 0.03, 0.03, 0.018, 0.017, 0.033, 0.023]
    beta_delta_t = [2.58e-3, 2.37e-3, 1.93e-3, 2.03e-3, 0, 1.12e-3, 2.75e-3, 1.76e-3, 2.61e-3, 1.7e-3, 2.02e-3, 9.94e-4,
                -4.87e-4, 0, 0, 0, 1.1e-3, 0, -1.4e-3, -1.78e-3, 0, -6.93e-4, 2.6e-4, -7.18e-4]
    delta_delta_t = [0.21, 0.39, 0.16, 0.289, 0.45, 0.33, 0.165, 0.36, 0, 0.26, 0, 0.19,
                0.13, 0.07, 0, 0.046, -0.15, 0.037, 0.22, 0.307, 0, 0.119, 0, 0.074]
    mu_delta_t = [-0.22, -0.59, -0.36, -0.45, -0.68, -0.59, -0.245, -0.48, -0.213, -0.52, -0.253, -0.42,
                0.0024, -0.102, 0.046, -0.083, 0.11, -0.066, -0.17, -0.66, -0.072, -0.05, -0.015, -0.005]
    
    # Coefficienstime and equation for tau (decay):
    
    alpha_tau = [0.047, 0.087, 0.054, 0.0403, 0, 0.035, 0.03, 0.057, 0.03, 0.0311, 0.05, 0.052,
                0.037, 0.0557, 0.029, 0.045, 0.029, 0.038, 0.031, 0.06, 0.04, 0.051, 0.024, 0.022]  
    beta_tau = [0, -1.89e-3, 5.37e-5, -1.26e-3, 0, -1.27e-3, 2.75e-3, -1.36e-3, 8.6e-4, -6.4e-4, 8.9e-4, 0,
                0, -8.2e-4, 8.0e-4, -5.46e-4, 0, -1.34e-3, 0, -1.45e-3, 9.4e-4, -1.41e-3, 0, -1.65e-3]
    delta_tau = [0.48, 0.58, 0.41, 0.387, 0.19, 0.19, 0.58, 0.63, 0.35, 0.44, 0.16, 0.12,
                0.39, 0.51, 0.25, 0.46, 0.36, 0.48, 0.34, 0.51, 0.25, 0.438, 0.303, 0.44]
    gamma_tau = [0.82, 0.58, 0.73, 0.58, 0, 0, 0, 0, 0, 0, 0, 0, 1.73, 1.63, 1.61, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    mu_tau = [-0.75, -0.87, -0.51, -0.372, -0.07, -0.03, -0.97, -0.96, -0.62, -0.55, -0.387, -0.166,
                -0.59, -0.68, -0.31, -0.55, -0.38, -0.39, -0.44, -0.60, -0.34, -0.368, -0.22, -0.19]
    avg_gamma = 0.15

    
    # Coefficienstime and equation for gamma (decay):
    alpha_gamma = [-0.032, -0.048, -0.044, -0.0403, -0.062, -0.061, -0.027, -0.024, -0.039, -0.037, -0.052, -0.066,
                -0.014, -0.015, -0.024, -0.031, -0.025, -2.67e-2, -0.0149, -0.0197, -0.028, -0.0334, -0.015, -0.0176] #<--should be =-0.048 for i=1? not =-0.48?
    beta_gamma = [-1.81e-3, -1.42e-3, -1.65e-3, -2.0e-3, -2.3e-3, -1.9e-3, -1.75e-3, -1.6e-3, -1.88e-3, -2.23e-3, -1.67e-3, -2.5e-3,
                -5.28e-4, -5.89e-4, -1.02e-3, -4.61e-4, -4.22e-4, 2.0e-4, -4.64e-4, 0, -8.32e-4, 0, 0, 5.65e-4]
    delta_gamma = [-0.1, -0.13, -0.16, 0, 0, 0.11, -0.18, -0.24, -0.18, -0.14, -0.21, 0,
                -0.11, -0.163, -0.055, -0.162, -0.145, -0.217, -0.122, -0.242, -0.123, -0.21, -0.229, -0.25]
    tau_gamma = [0.27, 0.26, 0.33, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0.38, 0.39, 0.36, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    mu_gamma = [0.64, 0.71, 0.72, 0.578, 0.61, 0.39, 0.74, 0.84, 0.76, 0.71, 0.849, 0.63,
                0.26, 0.299, 0.207, 0.302, 0.262, 0.274, 0.255, 0.378, 0.325, 0.325, 0.309, 0.236]
    avg_gamma = 0.15
    

    stat_err = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    sta_corr =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    # coefficienstime
    t_rise_p = 10**(alpha_t_rise[Pcoeff] * M + beta_t_rise[Pcoeff] * dist_in_km + delta_t_rise[Pcoeff] * log10(dist_in_km) + mu_t_rise[Pcoeff])
    t_rise_s = 10**(alpha_t_rise[Scoeff] * M + beta_t_rise[Scoeff] * dist_in_km + delta_t_rise[Scoeff] * log10(dist_in_km) + mu_t_rise[Scoeff])
    delta_t_p = 10**(alpha_delta_t[Pcoeff] * M + beta_delta_t[Pcoeff] * dist_in_km + delta_delta_t[Pcoeff] * log10(dist_in_km) + mu_delta_t[Pcoeff])
    delta_t_s = 10**(alpha_delta_t[Scoeff] * M + beta_delta_t[Scoeff] * dist_in_km + delta_delta_t[Scoeff] * log10(dist_in_km) + mu_delta_t[Scoeff])
    tau_p = 10**(alpha_tau[Pcoeff] * M + beta_tau[Pcoeff] * dist_in_km + delta_tau[Pcoeff] * log10(dist_in_km) + mu_tau[Pcoeff])
    tau_s = 10**(alpha_tau[Scoeff] * M + beta_tau[Scoeff] * dist_in_km + delta_tau[Scoeff] * log10(dist_in_km) + mu_tau[Scoeff])
    gamma_p = 10**(alpha_gamma[Pcoeff] * M + beta_gamma[Pcoeff] * dist_in_km + delta_gamma[Pcoeff] * log10(dist_in_km) + mu_gamma[Pcoeff])
    gamma_s = 10**(alpha_gamma[Scoeff] * M + beta_gamma[Scoeff] * dist_in_km + delta_gamma[Scoeff] * log10(dist_in_km) + mu_gamma[Scoeff])
    
    # Other variable (turn on saturation for larger evenstime?)
    C_p = (arctan(M-5) + (pi/2))*(c1[Pcoeff]*exp(c2[Pcoeff] * (M-5)))
    C_s = (arctan(M-5) + (pi/2))*(c1[Scoeff]*exp(c2[Scoeff] * (M-5)))
    R1 = sqrt(dist_in_km**2 + 9)
    
    # Basic AMplitudes
    A_p = 10**(a[Pcoeff]*M + b[Pcoeff]*(R1 + C_p) + d[Pcoeff]*log10(R1+C_p) + e[Pcoeff]+(sta_corr[Pcoeff]) + stat_err[Pcoeff])
    A_s = 10**(a[Scoeff]*M + b[Scoeff]*(R1 + C_s) + d[Scoeff]*log10(R1+C_s) + e[Scoeff]+(sta_corr[Scoeff]) + stat_err[Scoeff])
    
    # calculate envelope (ENV)
    envelope = zeros(len(times))

    # P envelope
    indx = where((times>=ptime) & (times<ptime+t_rise_p)) # between trigger and rise time
    if len(indx): envelope[indx] = (A_p/t_rise_p*(times[indx]-ptime)) # make sure we have data in that time frame and get envelope
    indx = where((times>=ptime+t_rise_p) & (times<ptime+t_rise_p+delta_t_p)) # flat area
    if len(indx): envelope[indx] = A_p # make sure we have data in that time frame and get envelope
    indx = where(times>ptime+t_rise_p+delta_t_p) # coda
    if len(indx): envelope[indx] = (A_p/((times[indx]-ptime-t_rise_p-delta_t_p+tau_p)**gamma_p)) # make sure we have data in that time frame and get envelope
    
    # S envelope
    indx = where((times>=stime) & (times<stime+t_rise_s)) # between trigger and rise time
    if len(indx): envelope[indx] += (A_s/t_rise_s*(times[indx]-stime)) # make sure we have data in that time frame and get envelope
    indx = where((times>=stime+t_rise_s) & (times<stime+t_rise_s+delta_t_s)) # flat area
    if len(indx): envelope[indx] += A_s # make sure we have data in that time frame and get envelope
    indx = where(times>stime+t_rise_s+delta_t_s) # coda
    if len(indx): envelope[indx] += (A_s/((times[indx]-stime-t_rise_s-delta_t_s+tau_s)**gamma_s)) # make sure we have data in that time frame and get envelope
    
    return envelope
    

def conically_avg_radiation_pattern(strike,dip,rake,azimuth,take_off_angle,
                                    component_angle,angle_range=45,Nrandom=1000):
    '''
    Get conically averaged radiation pattern, this is meant for horizontal 
    channels
    '''
    from numpy.random import rand
    from numpy import sin,cos,deg2rad,sign
    
    #Sample randomly over angle_range
    st=strike+angle_range*(rand(Nrandom)-0.5)
    di=dip+angle_range*(rand(Nrandom)-0.5)
    ra=rake+angle_range*(rand(Nrandom)-0.5)
    az=azimuth+angle_range*(rand(Nrandom)-0.5)
    to=take_off_angle+angle_range*(rand(Nrandom)-0.5)
 
    #Get theoretical radiation pattern 
    P_theor,SV_theor,SH_theor=radiation_pattern(strike,dip,rake,azimuth,take_off_angle)    
 
    #Get values of radiation pattern at specified angles
    P,SV,SH=radiation_pattern(st,di,ra,az,to)   
 
    #Projection angles to radians
    azimuth=deg2rad(azimuth)
    component_angle=deg2rad(component_angle)
    
    #Project SV and SH components onto component angle
    rad_pattern_theor = SV_theor*cos(component_angle-azimuth)+SH_theor*sin(component_angle-azimuth)
    rad_pattern = SV*cos(component_angle-azimuth)+SH*sin(component_angle-azimuth)
    
    #Square and sum
    rad_pattern=sum(rad_pattern**2)
    
    #Geometric mean
    rad_pattern=(rad_pattern/Nrandom)**0.5
    
    #What was polarity of theoretical rad. pattern?
    polarity=sign(rad_pattern_theor)
    
    #Apply
    rad_pattern=rad_pattern*polarity
    
    return rad_pattern
    

    
        
def conically_avg_P_radiation_pattern(strike,dip,rake,azimuth,take_off_angle,
                                    angle_range=45,Nrandom=1000):
    '''
    Get conically averaged radiation pattern, this is meant for horizontal 
    channels
    '''
    from numpy.random import rand
    from numpy import sin,cos,deg2rad,sign
    
    #Sample randomly over angle_range
    st=strike+angle_range*(rand(Nrandom)-0.5)
    di=dip+angle_range*(rand(Nrandom)-0.5)
    ra=rake+angle_range*(rand(Nrandom)-0.5)
    az=azimuth+angle_range*(rand(Nrandom)-0.5)
    to=take_off_angle+angle_range*(rand(Nrandom)-0.5)
 
    #Get theoretical radiation pattern 
    P_theor,SV_theor,SH_theor=radiation_pattern(strike,dip,rake,azimuth,take_off_angle)    
 
    #Get values of radiation pattern at specified angles
    P,SV,SH=radiation_pattern(st,di,ra,az,to)   
 
    #Projection angles to radians
    azimuth=deg2rad(azimuth)
    
    #Project SV and SH components onto component angle
    rad_pattern_theor = P_theor
    rad_pattern = P
    
    #Square and sum
    rad_pattern=sum(rad_pattern**2)
    
    #Geometric mean
    rad_pattern=(rad_pattern/Nrandom)**0.5
    
    #What was polarity of theoretical rad. pattern?
    polarity=sign(rad_pattern_theor)
    
    #Apply
    rad_pattern=rad_pattern*polarity
    
    return rad_pattern            
                
                    
                        
                            
                                    
def conically_avg_vert_radiation_pattern(strike,dip,rake,azimuth,take_off_angle,
                                    angle_range=40,Nrandom=1000):
    '''
    Get conically averaged radiation pattern, for the vertical channels
    '''
    from numpy.random import rand
    from numpy import sin,cos,deg2rad,sign,pi ,arccos,ones,rad2deg
  
    #get theoretical rad pattern
    P,SV,SH=radiation_pattern(strike,dip,rake,azimuth,take_off_angle) 
  
    #To radians
    rake=deg2rad(rake)                            
    strike=deg2rad(strike)
    dip=deg2rad(dip)
    azimuth=deg2rad(azimuth)
    take_off_angle=deg2rad(take_off_angle)                                                                                                                                  
    
    #Angles to use
    theta1=take_off_angle-deg2rad(angle_range)
    theta2=take_off_angle+deg2rad(angle_range)
      
    #Make sure angles aren't larger/smaller than reasonable
    if(theta1<pi/2):
        theta1=pi/2
    if(theta2>pi):
        theta2=pi 

    #Random numbers we'll use
    rand_num=rand(Nrandom)                                                                
                                               
    #Get suite of take off angles (theta) and azimuths (fa)                                                           
    theta=arccos((1.0-rand_num)*cos(theta1)+rand_num*cos(theta2))    
    fa=2*pi*rand(Nrandom) 
      
    #Make vector of angles
    strike=strike*ones(Nrandom)
    dip=dip*ones(Nrandom)
    rake=rake*ones(Nrandom)
           
    # Get radiation patterns                              
    P,SV,SH=radiation_pattern(rad2deg(strike),rad2deg(dip),rad2deg(rake),rad2deg(fa),rad2deg(theta))                           
    
    #Project
    SV=SV*sin(theta)
    
    #Sum
    SV=sum(abs(SV))
    
    #Get mean
    SVH=SV/(2*Nrandom)             
                                                  
    return SVH
    
    
    
def radiation_pattern(strike,dip,rake,azimuth,take_off_angle):
    '''
    Build theoretical P,SV and SH radiation patterns see Lay & Wallace p.340
    
    Inputs are in DEGREES
    
    '''
    
    from numpy import cos,sin,deg2rad
    
    #To radians
    rake=deg2rad(rake)                            
    strike=deg2rad(strike)
    dip=deg2rad(dip)
    azimuth=deg2rad(azimuth)
    take_off_angle=deg2rad(take_off_angle)    
                                                
    SR=sin(rake)                                                               
    CR=cos(rake)                                                               
    SD=sin(dip)                                                               
    CD=cos(dip)                                                               
    ST=sin(take_off_angle)                                                                
    CT=cos(take_off_angle)                                                                
    SS=sin(azimuth-strike)                                                            
    CS=cos(azimuth-strike)      
      
    # P wave                                                            
    P = CR*SD*ST**2*2*SS*CS - CR*CD*2*ST*CT*CS + SR*2*SD*CD*(CT**2-ST**2*SS**2) + SR*(CD**2-SD**2)*2*ST*CT*SS                                          
      
    #SV
    SV = SR*(CD**2-SD**2)*(CT**2-ST**2)*SS - CR*CD*(CT**2-ST**2)*CS + CR*SD*ST*CT*2*SS*CS - SR*SD*CD*2*ST*CT*(1+SS**2)  

    #SH                                    
    SH = CR*CD*CT*SS + CR*SD*ST*(CS**2-SS**2) + SR*(CD**2-SD**2)*CT*CS - SR*SD*CD*ST*2*SS*CS  

    return P,SV,SH
    
    
    
def get_P_wave_partition(incidence_angle,azimuth):
    '''
    Get partition of P-wave into horizontal and vertical components
    '''
    
    from numpy import sin,cos,deg2rad
    
    #Split into horizontal and vertical
    Zpartition=abs(cos(deg2rad(incidence_angle)))
    Hpartition=abs(sin(deg2rad(incidence_angle)))
    
    #Split horizontal into N and E
    Npartition=abs(Hpartition*cos(deg2rad(azimuth)))
    Epartition=abs(Hpartition*sin(deg2rad(azimuth)))
    
    return Npartition, Epartition, Zpartition
    
 
def is_subfault_in_smga(subfault_lon,subfault_lat,smga_center_lon,smga_center_lat,radius_in_km):
    '''
    Decide whether subfault is inside a defined strong motion generating area SMGA
    '''
    
    from pyproj import Geod
    
    p=Geod(ellps='WGS84')
    
    az,baz,d=p.inv(subfault_lon,subfault_lat,smga_center_lon,smga_center_lat)
    d=d/1000.
    
    if d<radius_in_km:
        subfault_in_smga=True
    else:
        subfault_in_smga=False
        
    return subfault_in_smga   