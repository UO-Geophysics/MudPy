def stochastic_simulation(home,project_name,rupture_name,GF_list,time_epi,model_name,
        rise_time_depths,hf_dt=0.01,stress_parameter=50e5,kappa=0.04,Qexp=0.6): 
    '''
    Run stochastic HF sims
    '''
    
    from numpy import genfromtxt,pi,logspace,log10,mean,where,exp
    from pyproj import Geod
    from obspy.geodetics import kilometer2degrees
    from obspy.taup import taup_create,TauPyModel
    from mudpy.forward import get_mu
    
    #Load the source
    fault=genfromtxt(rupture_name)    
    
    #Load stations
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    
    #load velocity structure
    structure=genfromtxt(home+project_name+'/structure/'+model_name)
    
    #Frequencies vector
    f=logspace(log10(hf_dt),log10(1/(2*hf_dt)),50)
    
    #Projection object for distance calculations
    g=Geod(ellps='WGS84')
    
    #Create taup velocity model object, paste on top of iaspei91
    taup_create.build_taup_model(home+project_name+'/structure/bbp_norcal.tvel',output_folder=home+project_name+'/structure/')
    velmod=TauPyModel(model=home+project_name+'/structure/bbp_norcal')
    taup_perturb=0.1
    
    #Moments
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    subfault_M0=slip*fault[:,10]*fault[:,11]*fault[:,13]
    M0=subfault_M0.sum()
    relative_subfault_M0=subfault_M0/M0
    
    #Corner frequency scaling
    i=where(slip>0)[0]
    N=len(i) #number of subfaults
    dl=mean(fault[:,10]**2+fault[:,11]**2)**0.5
    fc_scale=M0/(N*stress_parameter*dl**3)
    
    #Loop over stations
    for ksta in range(len(lonlat)):
    
        #Loop over subfaults
        for kfault in range(len(fault)):
            
            #Get subfault to station distance
            lon_source=fault[kfault,1]
            lat_source=fault[kfault,2]
            az,baz,dist=g.inv(lon_source,lat_source,lonlat[ksta,0],lonlat[ksta,1])
            dist_in_degs=kilometer2degrees(dist/1000.)
            
            #Get rho and beta at subfault depth
            zs=fault[kfault,3]
            mu,beta=get_mu(structure,zs,return_beta=True)
            rho=mu/beta**2
            
            #Get radiation scale factor (CURRENTLY MISSING CONICALLY AVERAGED RADIATION PATTERN)
            C=2./(4*pi*rho*beta**3)
            
            #Get local subfault rupture speed
            vr=get_local_rupture_speed(zs,beta,rise_time_depths)
            dip_factor=get_dip_factor(fault[kfault,5])
            
            #Subfault corner frequency
            fc_subfault=(2.1*vr)/(dip_factor*pi*dl)
            
            #get subfault source spectrum
            S=(relative_subfault_M0[kfault]*fc_scale*f**2)/(1+fc_scale*(f/fc_subfault)**2)
            
            #get high frequency decay
            P=exp(-pi*kappa*f)
            
            #get quarter wavelength amplificationf actors
            I=get_amplification_factors(f,structure,zs,beta,rho)
            
            #Get ray paths for all direct S arrivals
            paths=velmod.get_ray_paths(zs+taup_perturb,dist_in_degs,phase_list=['S','s'])
            
            #Get attenuation due to geometrical spreading (from the path length)
            path_length=get_path_length(paths[0],zs,dist_in_degs)
            
            #Get effect of intrinsic attenuation for that ray (path integrated)
            Q=get_attenuation(f,structure,paths[0],Qexp)
            
            #Build the entire path term
            G=(I/path_length)*Q
            
            #And finally multiply everything together to get the subfault amplitude spectrum
            A=C*S*G*P
            
                        
                        
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
    Get effect of intrinsic attenuation along the ray path
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
