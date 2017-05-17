def stochastic_simulation(home,project_name,rupture_name,GF_list,time_epi,model_name,
        rise_time_depths,total_duration=100,hf_dt=0.01,stress_parameter=50e5,
        kappa=0.04,Qexp=0.6,component='N'): 
    '''
    Run stochastic HF sims
    '''
    
    from numpy import genfromtxt,pi,logspace,log10,mean,where,exp,arange,zeros,argmin,ones
    from pyproj import Geod
    from obspy.geodetics import kilometer2degrees
    from obspy.taup import taup_create,TauPyModel
    from mudpy.forward import get_mu
    from obspy import Stream,Trace
    
    #initalize  output object
    st=Stream()
    
    #Load the source
    fault=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)    
    
    #Onset times for each subfault
    onset_times=fault[:,12]
    
    #Load stations
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[0],dtype='S')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    
    #load velocity structure
    structure=genfromtxt(home+project_name+'/structure/'+model_name)
    
    #Frequencies vector
    f=logspace(log10(hf_dt),log10(1/(2*hf_dt))+0.01,50)
    
    #Output time vector (0 is origin time)
    t=arange(0,total_duration,hf_dt)
    
    #Projection object for distance calculations
    g=Geod(ellps='WGS84')
    
    #Create taup velocity model object, paste on top of iaspei91
    taup_create.build_taup_model(home+project_name+'/structure/bbp_norcal.tvel',output_folder=home+project_name+'/structure/')
    velmod=TauPyModel(model=home+project_name+'/structure/bbp_norcal',verbose=True)
    taup_perturb=0.1
    
    #Moments
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    subfault_M0=slip*fault[:,10]*fault[:,11]*fault[:,13]
    subfault_M0=subfault_M0*1e7 #to dyne-cm
    M0=subfault_M0.sum()
    relative_subfault_M0=subfault_M0/M0
    
    #Corner frequency scaling
    i=where(slip>0)[0] #Non-zero faults
    N=len(i) #number of subfaults
    dl=mean((fault[:,10]+fault[:,11])/2) #perdominant length scale
    dl=dl/1000 # to km
    
    # Frankel 95 scaling of corner frequency
    fc_scale=(M0)/(N*50*dl**3*1e21) #Frankel scaling
    
    #Loop over stations
    for ksta in range(len(lonlat)):
    
        print '... working on '+component+' component semistochastic waveform for station '+sta[ksta]
    
        #initalize output seismogram
        tr=Trace()
        tr.stats.station=sta[ksta]
        tr.stats.delta=hf_dt
        tr.stats.starttime=time_epi
        hf=zeros(len(t))
        
        #Loop over subfaults
        for kfault in range(len(fault)):
            
            #Include only subfaults with non-zero slip
            if subfault_M0[kfault]>0:
            
                #Get subfault to station distance
                lon_source=fault[kfault,1]
                lat_source=fault[kfault,2]
                az,baz,dist=g.inv(lon_source,lat_source,lonlat[ksta,0],lonlat[ksta,1])
                dist_in_degs=kilometer2degrees(dist/1000.)
                
                #Get rho and beta at subfault depth
                zs=fault[kfault,3]
                mu,beta=get_mu(structure,zs,return_beta=True)
                rho=mu/beta**2
                
                #Get radiation scale factor
                if component=='N' or component=='E':
                    partition=1/2**0.5
                elif component=='Z':
                    partition=0.3 #this is addhoc right now, needs to be verified/adjsuted
                avg_radiation=0.63
                rho=rho/1000 #to g/cm**3
                beta=(beta/1000)*1e5 #to cm/s
                C=(2*partition*avg_radiation)/(4*pi*(rho)*(beta**3))
                
                #Get local subfault rupture speed
                beta=beta/100 #to m/s
                vr=get_local_rupture_speed(zs,beta,rise_time_depths)
                vr=vr/1000 #to km/s
                dip_factor=get_dip_factor(fault[kfault,5])
                
                #Subfault corner frequency
                fc_subfault=(2.1*vr)/(dip_factor*pi*dl)
                
                #get subfault source spectrum
                S=((relative_subfault_M0[kfault]*M0/N)*f**2)/(1+fc_scale*(f/fc_subfault)**2)
                
                #get high frequency decay
                P=exp(-pi*kappa*f)
                
                #get quarter wavelength amplificationf actors
                I=get_amplification_factors(f,structure,zs,beta,rho)
                
                #Get ray paths for all direct S arrivals
                Ppaths=velmod.get_ray_paths(zs-taup_perturb,dist_in_degs,phase_list=['P','p'])
                
                #Get ray paths for all direct S arrivals
                Spaths=velmod.get_ray_paths(zs+taup_perturb,dist_in_degs,phase_list=['S','s'])
                
                #Get attenuation due to geometrical spreading (from the path length)
                path_length=get_path_length(Spaths[0],zs,dist_in_degs)
                path_length=path_length/1000 #to km
                
                #Get effect of intrinsic aptimeenuation for that ray (path integrated)
                Q=get_attenuation(f,structure,Spaths[0],Qexp)
                
                #Build the entire path term
                G=(I/path_length)*Q
                
                #And finally multiply everything together to get the subfault amplitude spectrum
                A=C*S*G*P
                #A=S*G*P
                
                #Generate windowed time series
                duration=1./fc_subfault+0.063*(dist/1000)
                w=windowed_gaussian(duration,hf_dt,window_type='saragoni_hart')
                #w=windowed_gaussian(3*duration,hf_dt,window_type='cua',ptime=Ppaths[0].path['time'][-1],stime=Spaths[0].path['time'][-1])
                
                #Go to frequency domain, apply amplitude spectrum and ifft for final time series
                hf_seis=apply_spectrum(w,A,f,hf_dt)
                
                #What time does this time series start at?
                time_insert=Spaths[0].path['time'][-1]+onset_times[kfault]
                print 'ts = '+str(time_insert)+' , Td = '+str(duration)
                #time_insert=Ppaths[0].path['time'][-1]
                i=argmin(abs(t-time_insert))
                j=i+len(hf_seis)
                
                #Add seismogram
                hf[i:j]=hf[i:j]+hf_seis
                
        #Done add to trace and stream
        tr.data=hf/100 #convert to m/s**2
        st+=tr
    return st
            
                        
                        
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


def get_dip_factor(dip):
    if dip<45:
        dip_factor=0.82
    elif dip>60:
        dip_factor=1.0
    else:
        m=(0.18/15.)
        b=1.0-60*m
        dip_factor=m*dip+b
    return dip_factor                    
                                                                                                                        

                                                            
def get_amplification_factors(f,structure,zs,beta,rho):
    '''
    Get quarter wavelength amplificationf actors
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
    interpolator=interp1d(z,vs)
    interpolator2=interp1d(z,rho_model)
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
    interpolator=interp1d(fz,Afz)
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
    

def get_attenuation(f,structure,ray,Qexp):
    '''
    Get effect of intrinsic aptimeenuation along the ray path
    '''
    
    from numpy import diff,zeros,exp,pi
    from mudpy.forward import get_Q
    
    time=ray.path['time']
    time_in_layer=diff(time)
    depth=ray.path['depth']
    
    Qp=zeros(len(time_in_layer))
    Qs=zeros(len(time_in_layer))
    
    for k in range(len(Qp)):
        Qp[k],Qs[k]=get_Q(structure,depth[k])
        
    #Get the travel tiem weighted sum
    weightedQ=sum(time_in_layer/Qs)
    
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
        
    noise=noise*window
    
    return noise

        
def apply_spectrum(w,A,f,hf_dt):
    '''
    Apply the modeled spectrum to the windowed time series
    '''
    
    from numpy import fft,angle,cos,sin,sqrt,mean
    from scipy.interpolate import interp1d
    
    #to frequency domain
    fourier=fft.rfft(w)
    freq=fft.rfftfreq(len(w),hf_dt)
    
    #Make POWER spectrum of windowed time series have a mean of 1
    norm_factor=hf_dt*mean(abs(fourier)**2)**0.5
    #norm_factor=mean(abs(fourier)) 
    fourier=fourier/norm_factor
    
    #Keep phase
    phase=angle(fourier)
    
    #resample model amplitude spectr to frequencies
    interp=interp1d(f,A,bounds_error=False)
    amplitude=interp(freq)
    amplitude[0]=amplitude[1]/10
    
    #Apply model amplitude spectrum
    amplitude=amplitude*abs(fourier)
    
    #Obtain complex foureier series
    R=amplitude*cos(phase)
    I=amplitude*sin(phase)
    fourier=R+I*1j
    
    #ifft
    seis=fft.irfft(fourier)
    
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