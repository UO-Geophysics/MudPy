'''
D.Melgar 03/2016

Generate synthetic slip distributions using the K-L expansion method.
This grew (a lot) from RJ Leveque's KL2d_vonKarman notebook. I modified
it to work in a similar way to the MudPy ivnersions
'''


#Initalize project folders
def init(home,project_name):
    '''
    Initalizes file structure for a new problem
    
    IN:
        home: What dir will you be working from
        project_name: What name will you give this problem
        
    OUT:
        Nothing
    '''
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    clob='y'
    proj_dir=path.expanduser(home+project_name+'/')
    if path.exists(proj_dir):  #Path exists, clobber?
        clob=input('Project directory exists, clobber (y/n)?')
        if clob == 'y' or clob == 'Y': #Clobber baby
            clob=input('This will delete everything in this project directory, so, take a minute, think about it: clobber (y/n)?')
            if clob == 'y' or clob == 'Y':
                rmtree(proj_dir)
            else: #Leave direcory alone
                print('Phew, almost shot yourself in the foot there didn\'t you?')
        else: #Leave direcory alone
            print('Phew, almost shot yourself in the foot there didn\'t you?')
    if clob == 'y' or clob == 'Y':
        makedirs(proj_dir)
        #And make the subdirectories
        makedirs(proj_dir+'GFs')
        makedirs(proj_dir+'GFs/static')
        makedirs(proj_dir+'GFs/dynamic')
        makedirs(proj_dir+'GFs/matrices')
        makedirs(proj_dir+'GFs/STFs')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'data/distances')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'scripts')
        makedirs(proj_dir+'forward_models')
        makedirs(proj_dir+'output/ruptures')
        makedirs(proj_dir+'output/statics')
        makedirs(proj_dir+'output/waveforms')
        makedirs(proj_dir+'logs')
        makedirs(proj_dir+'analysis')
        makedirs(proj_dir+'analysis/frequency')
        #Copy templates into appropriate files
        mudpy=environ['MUD']+'/run/'
        #copy(mudpy+'template.fault',proj_dir+'data/model_info/')
        #copy(mudpy+'template.gflist',proj_dir+'data/station_info/')
        #copy(mudpy+'template.sta',proj_dir+'data/station_info/')
        #copy(mudpy+'template.mod',proj_dir+'structure/')


def llz2utm(lon,lat,projection_zone='None'):
    '''
    Convert lat,lon to UTM
    '''
    from numpy import zeros,where,chararray
    import utm
    from pyproj import Proj
    from scipy.stats import mode
    
    x=zeros(lon.shape)
    y=zeros(lon.shape)
    zone=zeros(lon.shape)
    b=chararray(lon.shape)
    if projection_zone==None:
        #Determine most suitable UTM zone
        for k in range(len(lon)):
            #x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k]-360)
            x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k])
        zone_mode=mode(zone)
        i=where(zone==zone_mode)[0]
        letter=b[i[0]]
        z=str(int(zone[0]))+letter
    else:
        z=projection_zone
    print(z)
    p = Proj(proj='utm',zone=z,ellps='WGS84')
    x,y=p(lon,lat)
    return x,y



def subfault_distances_3D(home,project_name,fault_name,slab_name,projection_zone):
    """
    Estimate the distance between subfaults i and j for every pair in the list
    fault.subfaults. For a 3D fault geometry

    :Inputs:
      -  *fault* of MudPy .fault class
    
    :Outputs:
      - *D* array of Euclidean distances based on longitudes, latitudes, and depths
      - *Dstrike* array of estimated distances along strike direction
      - *Ddip* array of estimated distances along dip direction
    with D**2 = Dstrike**2 + Ddip**2 to within roundoff.

    For each array, the [i,j] entry is distance from subfault i to j when
    ordered in the order the subfaults appear in the list fault.subfaults.

    Distance in dip direction based on differences in depth and average dip of the fault.  

    """

    from numpy import sqrt,sin,cos,deg2rad,zeros,meshgrid,linspace,where,c_,unravel_index,sort,diff,genfromtxt,sign,argmin
    from scipy.interpolate import griddata
    from matplotlib import pyplot as plt
    from scipy.spatial.distance import cdist
    from pyproj import Geod
    
    #if you want the simplified distances
    if slab_name==None:
        #Read fault geometry data
        fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
        
        #Initalize distance output arrays
        nsubfaults = len(fault)
        Dstrike = zeros((nsubfaults,nsubfaults))
        Ddip = zeros((nsubfaults,nsubfaults))  
        
        #What's the average dip of the fault model?? 
        dip=fault[:,5].mean()
        
        #What's the average strike of the fault model?? 
        strike=fault[:,4].mean()
        
        #Instantiate projection object
        g = Geod(ellps='WGS84') 
        
        #Loop over faults an compute distances
        print('Getting inter-fault distances')
        for i in range(len(fault)):
            
            if i%10==0:
                print('... working on subfault '+str(i)+' of '+str(len(fault)))
            
            #For each subfault loop over every other subfault
            for j in range(len(fault)):
                
                #Subfault distance with itslef is zero
                if i==j:
                    Ddip[i,j] = 0
                    Dstrike[i,j] = 0
                else:
                    lon_origin=fault[i,1]
                    lat_origin=fault[i,2]
                    lon_target=fault[j,1]
                    lat_target=fault[j,2]
                    az,baz,dist=g.inv(lon_origin,lat_origin,lon_target,lat_target)
                    
                    #straight line distance
                    delta_strike=dist/1000
                    #define projection angel as azimuth minus average strike
                    alpha=az-strike
                    #Project distance
                    delta_strike=abs(delta_strike*cos(deg2rad(alpha)))
                    
                    #Down dip is jsut depth difference / avg dip of model
                    z_origin=fault[i,3]
                    z_target=fault[j,3]
                    delta_dip=abs(z_origin-z_target)/sin(deg2rad(dip))
                    
                    #get the correct signs
                    if z_origin<z_target:
                        dip_sign=-1.
                    else:
                        dip_sign=1.
                        
                        
                    #Fix azimuth polarity
                    if az<0:
                        az+=360
                        
                    #Rotate and rectify azimuth
                    azimuth_rotated=az-strike
                    if azimuth_rotated<0:
                        azimuth_rotated+=360
                    
                        
                    if azimuth_rotated >=90 and azimuth_rotated <= 270:
                        strike_sign = -1.
                    else:
                        strike_sign= 1.
                    
                    Ddip[i,j] = dip_sign*delta_dip
                    Dstrike[i,j] = strike_sign*delta_strike
           
                    
    #If there's a slab_model file and you want the onfault complicated to get distances
    else:
        
        #Load things
        fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
        slab_model=genfromtxt(home+project_name+'/data/model_info/'+slab_name)    
    
        #Initalize distance output arrays
        nsubfaults = len(fault)
        Dstrike = zeros((nsubfaults,nsubfaults))
        Ddip = zeros((nsubfaults,nsubfaults))
        
        #Get average dip
        avg_dip=fault[:,5].mean()
        
        #get down-dip azimuths
        down_dip=fault[:,4]+90
        # i=where(down_dip>360)[0]
        # down_dip[i]=down_dip[i]-360
        
        #Convert slab1.0 to local UTM coordinates
        slab_x,slab_y=llz2utm(slab_model[:,0],slab_model[:,1],projection_zone)
        slab_x,slab_y = slab_x/1000,slab_y/1000
        slab_z=-slab_model[:,2]
        
        #Convert faul centroid coordinates to local UTM
        fault_x,fault_y=llz2utm(fault[:,1],fault[:,2],projection_zone)
        fault_x,fault_y = fault_x/1000,fault_y/1000
        fault_z=fault[:,3]
        
        #The goal is to only keep slab points close tot he actual fault, to accomplish this
        # Calcualte the distance from each fault to slab point and only keep things
        #Within a minimum distance
        slab2fault_min_distance=30 #in km
        keep_slab=zeros((len(fault_x),len(slab_x)))
        for ksub in range(len(fault)):
            dist=((fault_x[ksub]-slab_x)**2+(fault_y[ksub]-slab_y)**2)**0.5
            i=where(dist<slab2fault_min_distance)[0]
            keep_slab[ksub,i]=1
        i=where(keep_slab.sum(axis=0)>0)[0]
        slab_x=slab_x[i]
        slab_y=slab_y[i]
        slab_z=slab_z[i]
        
        
        # grid Slab1.0 for making depth contours to be then used for along-strike distance calculation
        ngrid_pts=500
        X=linspace(slab_x.min(),slab_x.max(),ngrid_pts)
        Y=linspace(slab_y.min(),slab_y.max(),ngrid_pts)
        X,Y = meshgrid(X,Y)
        Z = griddata(c_[slab_x, slab_y], slab_z, (X, Y), method='linear')
        
        # x-coordinates for down_dip line
        x_range=slab_x.max()-slab_x.min()
        x_down_dip=linspace(-x_range/2,x_range/2,200)
        
        #get contours
        all_contours=[]
        contour_lengths=zeros(len(fault))
        print('Calculcating slab contours')
        for i in range(len(fault)):
            if i%10==0:
                print('... working on subfault '+str(i)+' of '+str(len(fault)))
            
            contour=plt.contour(X,Y,Z,levels=[fault[i,3]])
            contour=contour.collections[0].get_paths()[0].vertices
            contour_lengths[i]=((contour[0,0]-contour[-1,0])**2+(contour[0,1]-contour[-1,1])**2)**0.5
            all_contours.append(contour)
        
        #if a contour is shorter than this number use the next deepest
        minimum_contour_length=0.95*contour_lengths.max()
        
        
        #Loop over number of subfaults, we want the distance from i-th fault to all other (j) subfaults
        print('Getting inter-fault distances')
        for i in range(len(fault)):
            if i%10==0:
                print('... working on subfault '+str(i)+' of '+str(len(fault)))
            #Current fault
            xi = fault_x[i]
            yi = fault_y[i]
            zi = fault_z[i]
            
            #find most approriate contour for lenght calculation
            icontour=where(contour_lengths>minimum_contour_length)[0]
            
            #find closest depth_contour of the ones that pass the minimum length
            deltaZ=abs(fault_z[i]-fault_z)
            icontour_depth=argmin(deltaZ[icontour])
            
            #This si the contour that is long enough and closest in depth
            icontour_correct=icontour[icontour_depth]
            contour=all_contours[icontour_correct]
            
            
#            # X, Y and Z are matrices with the grid info, now create one contour at each subfault centroid depth
#            contour=plt.contour(X,Y,Z,levels=fault[i,3])
#            
#            #Get contour at depth of current subfault
#            contour=contour.collections[0].get_paths()[0].vertices
            
            # Now find coordinates of point on this contour closest to subfault centroid
            dist=sqrt((xi-contour[:,0])**2+(yi-contour[:,1])**2)
            imin=dist.argmin()
            
            # These are coordinates on the contour
            xmin_i=contour[imin,0]
            ymin_i=contour[imin,1]
            
            #For each subfault loop over every other subfault
            for j in range(len(fault)):
                xj = fault_x[j]
                yj = fault_y[j]
                zj = fault_z[j]
                
                #Get down_dip y coordinates
                y_down_dip=x_down_dip*(-cos(deg2rad(down_dip[j])))
                
                #Move line origin to subfault being tested
                x_down_dip_subfault=x_down_dip+xj
                y_down_dip_subfault=y_down_dip+yj
                
                #Get coordinates of intersection point between contour and down-dip line by finding minimum distance
                dist=cdist(contour,c_[x_down_dip_subfault,y_down_dip_subfault])
                r,c=unravel_index(dist.argmin(),dist.shape)
                xmin_j=contour[r,0]
                ymin_j=contour[r,1]           
                
                #Keep only points on the contour array that correspond to starting and stopping points along the path
                keep=sort([r,imin])
                contour_integral=contour[keep[0]:keep[1]+1]
                
                #Along strike distance is the path integral along contour between (xmin_i,ymin_i) and (xmin_j,ymin_j)
                dx=diff(contour_integral[:,0])
                dy=diff(contour_integral[:,1])
                delta_strike=sqrt(dx**2+dy**2).sum() #Phew, that was hard
                #Give negative sign if subfault is down_strike
                strike_sign=sign(ymin_j-ymin_i)
                delta_strike=strike_sign*delta_strike
                
    
                
                #get down dip distance from depth and average dip
                delta_dip=(zi-zj)/sin(deg2rad(avg_dip))
            
                #Now the outputs
                if i==j:
                    Ddip[i,j] = 0
                    Dstrike[i,j] = 0
                else:
                    Ddip[i,j] = delta_dip
                    Dstrike[i,j] = delta_strike
                
    return Dstrike,Ddip



def get_mean_slip(target_Mw,fault_array,vel_mod):
    '''
    Depending on the target magnitude calculate the necessary uniform slip
    on the fault given the 1D layered Earth velocity model
    '''
    
    from numpy import genfromtxt,zeros,ones
    from mudpy.forward import get_mu
    
    vel=genfromtxt(vel_mod)
    areas=fault_array[:,8]*fault_array[:,9]
    mu=zeros(len(fault_array))
    for k in range(len(mu)):
        mu[k]=get_mu(vel,fault_array[k,3])
    target_moment=10**(1.5*target_Mw+9.1)
    mean_slip=ones(len(fault_array))*target_moment/sum(areas*mu)
    
    return mean_slip,mu
    

    
def vonKarman_correlation(Dstrike,Ddip,Lstrike,Ldip,hurst):
    '''
    Create VonKarman correlation function
    
    Lstrike,Ldip are correlation lengths
    hurst is the Hurst exponent
    '''
    from scipy.special import kv as Bessel # Modified Bessel function of second kind
    from numpy import sqrt
    
    delta = 1e-7  # fudge factor to avoid divide by zero
    
    #Get rid of negatives
    Dstrike=abs(Dstrike)
    Ddip=abs(Ddip)
    
    #Define von Karmen ACF in terms of G_H(r) where H is the Hurst exponent:
    G = lambda r,H: (r+delta)**H * Bessel(H,(r+delta))
    vonKarman = lambda r,H: G(r,H)/G(delta,H)
    
    r = sqrt((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2)
    C = vonKarman(r,hurst)
    
    return C


def get_lognormal(mean_slip,C,target_Mw,fault_array,vel_mod,alpha=0.6):
    '''
    Exponentiate to get lognormal correlation
    '''
    
    from numpy import log,diag

    #Get sigma as a fraction of mean slip
    sigma_slip=alpha*mean_slip
    #Now generate desired correlation
    Cov_g = log((sigma_slip/mean_slip) * (C*(sigma_slip/mean_slip)).T + 1.)
    mean_slip_g = log(mean_slip) - diag(Cov_g)/2.
    
    return Cov_g,mean_slip_g

     
def get_covariance(mean_slip,C,target_Mw,fault_array,vel_mod,alpha=0.6):
    '''
    Exponentiate to get lognormal correlation
    '''
    
    #Get sigma as a fraction of mean slip
    sigma_slip=alpha*mean_slip
    #Now generate desired covariance
    Cov = sigma_slip * (C*sigma_slip).T
    
    return Cov
           
                             
def get_eigen(C):
    '''
    Get eigen vectors/values from correlation matrix and return sorted results
    '''
    from numpy.linalg import eig
    from numpy import real,argsort
    
    eigenvals, V = eig(C)
        
    eigenvals = real(eigenvals)  # imaginary parts should be at rounding level
    V = real(V)
    
    # Sort eigenvalues:
    i = list(argsort(eigenvals))
    i.reverse()
    eigenvals = eigenvals[i]
    V = V[:,i]
    
    return eigenvals,V
    
    
def make_KL_slip(fault,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=True,maxiter=5,seed=12345):
    '''
    Make slip map using num_modes
    '''
    from numpy import sqrt,exp
    from numpy.random import randn
    from numpy.random import seed as random_seed
    
    iterations=0
    success=False
    

    
    while True:
        #Generate random numbers
            #Is there a seed?
        if seed != None:
            random_seed(seed)
        if len(fault)>num_modes:
            z = randn(num_modes) 
        else: #if fewer faults than requested modes then use all modes
            z = randn(len(fault))
        KL_slip = mean_slip.copy()  # start with the mean slip
        # add in the terms in the K-L expansion:
        for k in range(len(z)):
            if eigenvals[k]>0: #only use positive eigenvals
                KL_slip += z[k] * sqrt(eigenvals[k]) * V[:,k]
        # exponentiate for lognormal:
        if lognormal==True:
            KL_slip = exp(KL_slip)
        #Check if max_slip condition is met, if so then you're done
        if KL_slip.max()<=max_slip:
            success=True
            break
        iterations+=1
        if iterations>maxiter:
            print('... ... ... improper eigenvalues, recalculating...')
            break
    

    
    return KL_slip,success


def rectify_slip(slip_unrectified,percent_reject=10):
    '''
    Deal with negative slip values
    '''
    
    from numpy import where
    
    slip=slip_unrectified.copy()
    #Find the negatives
    i=where(slip_unrectified<0)[0]
    percent_negative=(float(len(i))/len(slip_unrectified))*100
    if percent_negative>percent_reject: #Too many negatives
        rejected=True
    else:
        slip[i]=0
        rejected=False
    
    return slip,rejected,percent_negative



def select_faults(whole_fault,Dstrike,Ddip,target_Mw,num_modes,scaling_law,
    force_area,no_shallow_epi=True,hypo_depth=10,param_norm=(0.0451,0.1681),no_random=False,
    subfault_hypocenter=None,use_hypo_fraction=True,option=0):
    '''
    Select a random fault to be the hypocenter then based on scaling laws and a 
    target magnitude select only faults within the expected area plus some 
    buffer factor
    '''
    
    from numpy.random import randint,normal
    from numpy import array,where,argmin,arange,log10
    from scipy.stats import norm,expon
    from random import choice
    
    
    #Select random subfault as center of the locus of L and W
    if subfault_hypocenter is None: 
        hypo_fault=randint(0,len(whole_fault)-1)
    else: #get subfault closest to hypo
        hypo_fault=subfault_hypocenter
    
    if force_area==True and no_random==False: #Use entire fault model  nothing more to do here folks
        
        selected_faults =  arange(len(whole_fault))
        Lmax=Dstrike[selected_faults,:][:,selected_faults].max()    
        Wmax=Ddip[selected_faults,:][:,selected_faults].max()
        #Convert to effective length/width
        Leff=0.85*Lmax
        Weff=0.85*Wmax
        
        return selected_faults,hypo_fault,Lmax,Wmax,Leff,Weff,0,0,0 #extra zeros after Yu-Sheng's updates, see return call below
        
    
    if force_area==False and no_random==True: #Use the Blasser medain L and W
        if scaling_law.upper()=='T':
            length=10**(-2.37+0.57*target_Mw)
            width=10**(-1.86+0.46*target_Mw)
        elif scaling_law.upper()=='S':
            length=10**(-2.69+0.64*target_Mw)
            width=10**(-1.12+0.3*target_Mw) 
        elif scaling_law.upper()=='N':
            length=10**(-1.91+0.52*target_Mw)
            width=10**(-1.2+0.36*target_Mw)
    
    if force_area==False and no_random==False: #Use scaling laws from Blaser et al 2010
        if scaling_law.upper()=='T':
            length_mean=-2.37+0.57*target_Mw
            length_std=0.18
            length=10**normal(length_mean,length_std)
            width_mean=-1.86+0.46*target_Mw
            width_std=0.17
            width=10**normal(width_mean,width_std)
        elif scaling_law.upper()=='S':
            length_mean=-2.69+0.64*target_Mw
            length_std=0.18
            length=10**normal(length_mean,length_std)
            width_mean=-1.12+0.3*target_Mw
            width_std=0.15
            width=10**normal(width_mean,width_std) 
        elif scaling_law.upper()=='N':
            length_mean=-1.91+0.52*target_Mw
            length_std=0.18
            length=10**normal(length_mean,length_std)
            width_mean=-1.2+0.36*target_Mw
            width_std=0.16
            width=10**normal(width_mean,width_std)           
        elif scaling_law.upper()=='SSE':
            if target_Mw>6.74: # Mw=6.74
                length_mean=-4.4848+1.0564*target_Mw
                length_std=0.1506
                length=10**normal(length_mean,length_std)
                width_mean=-0.4532+0.3593*target_Mw
                width_std=0.0930
                width=10**normal(width_mean,width_std) # Schmidt & Gai 2010, M - T
                option=1
            elif target_Mw<5.88: # Mw=5.80
                length_mean=-2.5864+0.7948*target_Mw
                length_std=0.0849
                length=10**normal(length_mean,length_std)
                width_mean=-0.2750+0.3163*target_Mw
                width_std=0.0653
                width=10**normal(width_mean,width_std) # Michel et al. 2019, M - T^3
                option=2
            else:
                option=[1,2]
                option=choice(option)
                if option==1:
                    length_mean=-4.4848+1.0564*target_Mw
                    length_std=0.1506
                    length=10**normal(length_mean,length_std)
                    width_mean=-0.4532+0.3593*target_Mw
                    width_std=0.0930
                    width=10**normal(width_mean,width_std) # Schmidt & Gai 2010
                elif option==2:
                    length_mean=-2.5864+0.7948*target_Mw
                    length_std=0.0849
                    length=10**normal(length_mean,length_std)
                    width_mean=-0.2750+0.3163*target_Mw
                    width_std=0.0653
                    width=10**normal(width_mean,width_std) # Michel et al. 2019
            if width>60:
                  option1=[3,4]
                  option1=choice(option1)
                  option=option*option1
                  if option1==4:
                      a=length*width
                      width=60+10**normal(0,width_std)
                      length=a/width
        
    # #so which subfault ended up being the middle?
    # center_subfault=hypo_fault  # I'm not sure why this is getting defined here?
        
    # hypo_fault=int(hypo_fault)
    
    #Get max/min distances from hypocenter to all faults
    dstrike_max=Dstrike[:,hypo_fault].max()
    dstrike_min=Dstrike[:,hypo_fault].min()
    ddip_max=Ddip[:,hypo_fault].max()
    ddip_min=Ddip[:,hypo_fault].min()
    
    #Work on strike first
#    strike_bounds=array([0,length/2])
    
    strike_bounds=array([-length/2,length/2])
    
    if strike_bounds[0]<dstrike_min:#Length is outside domain
        strike_bounds[1]=strike_bounds[1]+abs(dstrike_min-strike_bounds[0])
        strike_bounds[0]=dstrike_min
    if strike_bounds[1]>dstrike_max:#Length is outside domain
        strike_bounds[0]=strike_bounds[0]-abs(dstrike_max-strike_bounds[1])
        strike_bounds[1]=dstrike_max
        
    #Now get dip ranges
#    dip_bounds=array([0,width/2])
    dip_bounds=array([-width/2,width/2])
    
    if dip_bounds[0]<ddip_min:#Length is outside domain
        dip_bounds[1]=dip_bounds[1]+abs(ddip_min-dip_bounds[0])
        dip_bounds[0]=ddip_min
    if dip_bounds[1]>ddip_max:#Length is outside domain
        dip_bounds[0]=dip_bounds[0]-abs(ddip_max-dip_bounds[1])
        dip_bounds[1]=ddip_max
    
    Ds=Dstrike[:,hypo_fault]
    Dd=Ddip[:,hypo_fault]
    
    #Now select faults within those distances
    selected_faults=where((Ds>=strike_bounds[0]) & (Ds<=strike_bounds[1]) & (Dd>=dip_bounds[0]) & (Dd<=dip_bounds[1]))[0]
    
    ## This code block is for selecting any random subfault as the hypocenter as oppsoed to something
    ## based on lieklihood of strike fraction and dip fraction
    ##From within the selected faults randomly select the hypocenter
    hypo_found=False
    

    if len(selected_faults)>1: #there is more than 1 subfault
    
        i=randint(0,len(selected_faults)-1)
        hypo_fault=selected_faults[i]
        if no_shallow_epi==True:
            while hypo_found==False:
                if whole_fault[hypo_fault,3]<hypo_depth:
                    print('... ... ... hypocenter is km too shallow at %dkm, recalculating...' %(whole_fault[hypo_fault,3]))
                    i=randint(0,len(selected_faults)-1)
                    hypo_fault=selected_faults[i]
                else:
                    hypo_found=True
                    
    else: #only one subfault
        hypo_fault=selected_faults[0]
        hypo_found=True
    ############
    
    
    
    
    # This code block is for selecting the hypocenter based on probability
    # of strike fraction and dip fraction, randomly generate it and find closest 
    # Ds and Dd
    
    if use_hypo_fraction == True:
        
        #Need to find "center" of the selected faults. To do this look for the subfault
        #with the lowest combined maximum along strike and along dip distance to all
        #other faults
        
        
        Dsmin=Dstrike[selected_faults,:] ; Dsmin=Dsmin[:,selected_faults]
        Ddmin=Ddip[selected_faults,:] ; Ddmin=Ddmin[:,selected_faults]
        
        min_dist=abs(Dsmin).max(axis=1) + abs(Ddmin).max(axis=1)
        imin=argmin(min_dist)
        center_subfault=selected_faults[imin]
        
        #get distances from everyone to "center subfault"
        Ds=Dstrike[center_subfault,selected_faults]
        Dd=Ddip[selected_faults,center_subfault]
        
        
        #Along dip parameters (normal distr)
        mu= -0.117
        sigma=0.18
        #along strike parameters (expon distr)
        lamb = 6.211
        
        
        
        while True:
            dip_fraction=norm.rvs(mu,sigma)
            if dip_fraction<0.5 and dip_fraction>-0.5:
                break
        
        while True:
            strike_fraction=expon.rvs(0,1/lamb)
            if strike_fraction<0.5:
                break
            
        # go from fractions to distances, *2 is necessary to use the whole 
        # distance range since the PDF int he Melgar&Hayes paper only goes to 
        # a fraction value of 0.5
        if dip_fraction<0:
            dip_distance=Dd.min()*(abs(dip_fraction)*2)
        else:
            dip_distance=Dd.max()*(dip_fraction*2)
            
        strike_distance=abs(Ds).max()*strike_fraction*2
        sign=2*randint(0,2)-1
        strike_distance*=sign
            
        
        
        #where is the fualt that is this distance from the middle?
        hypo_fault=argmin((Ds-strike_distance)**2+(Dd-dip_distance)**2)      
        hypo_fault=selected_faults[hypo_fault]

    #############
                
    #From the selected faults determine the actual along strike length (Leff) and down-dip width (Weff)
    #Check it doesn't exceed physically permissible thresholds
    Lmax=Dstrike[selected_faults,:][:,selected_faults].max()    
    Wmax=Ddip[selected_faults,:][:,selected_faults].max()
    
    #Convert to effective length/width
    Leff=0.85*Lmax
    Weff=0.85*Wmax
    
    
    return selected_faults,hypo_fault,Lmax,Wmax,Leff,Weff,option,length,width

    
        
def get_rise_times(M0,slip,fault_array,rise_time_depths,stoc_rake,rise_time='MH2017',rise_time_std=0.1,option=0):
    '''
    Calculate individual subfault rise times
    '''     
    
    from numpy import diff,ones,where,exp
    from numpy.random import randn,normal
    from numpy import arctan2,rad2deg,zeros
    from random import choice
    
    #Moment to dyne-cm (Because old seismologists...)
    M0=M0*1e7
    
    #Determine average rise time based on total moment of the event (Graves,Pitarka, 2010, eq. 8)
    if rise_time=='GP2010':
        tau_average=1.6*1e-9*M0**(1./3) #This is what Graves and Pitarka use in GP 2010
    elif rise_time=='GP2015':
        tau_average=1.45*1e-9*M0**(1./3) #This is GP2015
    elif rise_time=='S1999':
        tau_average=2.0*1e-9*M0**(1./3)  #This is the original from Sommerville 1999 SRL, page 74
    elif rise_time=='MH2017':
        tau_average=4.226*1e-8*M0**(0.293) #This is  Melgar & Hayes 2017
    elif rise_time=='SSE': # pending
    #tau_average=0.0145*M0**(1./3)
    # print(option)
    # tau_average=12.4792*M0**(0.1959)
    # rise_time_std=0.1454
    # tau_average=0.0039*M0**(1./3)
    # rise_time_std=0.1510
    
    # tau_average=0.0039*slip**(1./3)
    # rise_time_std=0.1510
        tau_average=0.0145*M0**(1./3)*0.4939#0.5984 # Michel et al. 2019
        rise_time_std=0.1429
    # if option==2:#M0<10**24.8: # Mw=5.80
    #     tau_average=0.0145*M0**(1./3)*0.4939#0.5984 # Michel et al. 2019
    #     rise_time_std=0.1429
    # else:#M0>10**26.2: # Mw=6.74
    #     tau_average=7.6155*1e-20*M0*0.4939#0.5984 # Schmidt & Gai 2010
    #     rise_time_std=0.1615
    # print(option)    
    
    #Determine slope and intercept of k-scaling line
    slope=1./(rise_time_depths[0]-rise_time_depths[1])
    intercept=1-slope*rise_time_depths[1]
    
    #Get rakes
    rake=stoc_rake.copy()
    i=where(rake<0)[0]
    rake[i]=360+rake[i]
    rake=rake.mean()
    
    #Apply GP 2015 FR value
    if (rake>0 and rake<180):
        FR=1-(rake-90)/90.
    else:
        FR=0

    #GP 2015 FD value
    dip=fault_array[:,5]
    dip=dip.mean()
    if dip >45:
        FD=1-(dip-45)/45.
    else:
        FD=1
    
    #GP 2015 alpha_t 
    alpha=1./(1+FD*FR*0.1)
    
    #rescale average rise time
    tau_average=tau_average*alpha
        
    #For each depth determine the value of depth scaling (this is for GP2010 eq.7)
    depth_scale=ones(len(fault_array))
    ishallow=where(fault_array[:,3]<=rise_time_depths[0])[0]
    depth_scale[ishallow]=2
    itransition=where((fault_array[:,3]>rise_time_depths[0]) & (fault_array[:,3]<rise_time_depths[1]))[0]
    depth_scale[itransition]=slope*fault_array[itransition,3]+intercept
    
    #Now determine the scaling constant k
    k=(len(slip)*tau_average)/(sum(depth_scale*slip**0.5))
    
    #Stochastic perturbations
    rand_num=randn(len(slip))
    perturbations=exp(rise_time_std*rand_num)

    #And on to the actual subfault rise times
    rise_times=perturbations*depth_scale*k*(slip**0.5)
    #rise_times=depth_scale*k*(slip**0.5)
    
    return rise_times
    
 

    
    
def get_rupture_onset(home,project_name,slip,fault_array,model_name,hypocenter,shypo,
        rise_time_depths,M0,velmod,sigma_rise_time=0.2,shear_wave_fraction_shallow=0.49,shear_wave_fraction_deep=0.8):
    
    home,project_name,slip,fault_array,model_name,hypocenter,rise_time_depths,M0,velmod,shear_wave_fraction_shallow,shear_wave_fraction_deep
    '''
    Using a custom built tvel file ray trace from hypocenter to determine rupture
    onset times
    '''
        
    from numpy import genfromtxt,zeros,arctan2,sin,r_,where,log10,isnan,argmin,setxor1d,exp
    from numpy .random import rand,randn,randint
    from obspy.geodetics import gps2dist_azimuth, kilometer2degrees
    import warnings
    import numpy as np
    
    # swf_r=randint(1,11)
    # shear_wave_fraction_shallow=1/60/60/24*swf_r
    # shear_wave_fraction_deep=1/60/60/24*swf_r
    
    #I don't condone it but this cleans up the warnings
    warnings.filterwarnings("ignore")
    
    #Load velocity model
    vel=genfromtxt(home+project_name+'/structure/'+model_name)
        
    # Convert from thickness to depth to bottom of layer
    depth_to_top=r_[0,vel[:,0].cumsum()[0:-1]]
        
    #Get rupture speed shear-wave multipliers
    rupture_multiplier=zeros(len(vel))
    # Shallow 
    i=where(depth_to_top<=rise_time_depths[0])[0]
    rupture_multiplier[i]=shear_wave_fraction_shallow
    # Deep 
    i=where(depth_to_top>=rise_time_depths[1])[0]
    rupture_multiplier[i]=shear_wave_fraction_deep
    # Transition 
    i=where((depth_to_top<rise_time_depths[1]) & (depth_to_top>rise_time_depths[0]))[0]
    slope=(shear_wave_fraction_deep-shear_wave_fraction_shallow)/(rise_time_depths[1]-rise_time_depths[0])
    intercept=shear_wave_fraction_deep-slope*rise_time_depths[1]
    rupture_multiplier[i]=slope*depth_to_top[i]+intercept
    
    
    #Perturb depths of the hypocenter so that faults at the same depth are not zero onset
    delta=0.00001
    i_same_as_hypo=where(fault_array[:,3]==hypocenter[2])[0]
    dist=((fault_array[:,1]-hypocenter[0])**2+(fault_array[:,2]-hypocenter[1])**2)**0.5
    i_hypo=argmin(dist)

    #Get faults at same depth that are NOT the hypo
    i_same_as_hypo=setxor1d(i_same_as_hypo,i_hypo)
    #perturb
    R=rand(1)[0]
    fault_array[i_hypo,3]=fault_array[i_hypo,3]-delta*R
    hypocenter[2]=hypocenter[2]-delta*R
    R=rand(len(i_same_as_hypo))
    fault_array[i_same_as_hypo,3]=fault_array[i_same_as_hypo,3]+delta*R
    
    #Loop over all faults
    t_onset=zeros(len(slip))
    length2fault=zeros(len(slip))
    #Perturb all subfault depths a tiny amount by some random number so that they NEVER lie on a layer interface
    z_perturb=(rand(len(fault_array))-0.5)*1e-6
    fault_array[:,3]=fault_array[:,3]+z_perturb
    for kfault in range(len(slip)):
        D,az,baz=gps2dist_azimuth(hypocenter[1],hypocenter[0],fault_array[kfault,2],fault_array[kfault,1])
        D=D/1000
        #Start and stop depths
        if fault_array[kfault,3]<=hypocenter[2]:
            zshallow=fault_array[kfault,3]
            zdeep=hypocenter[2]
        else:
            zdeep=fault_array[kfault,3]
            zshallow=hypocenter[2]
        #Get angle between depths
        theta=arctan2(zdeep-zshallow,D)
        # get hypotenuse distance on all layers
        delta_ray=vel[:,0]/sin(theta)
        # Calculate distance in each layer
        depth1=0
        depth2=vel[0,0]
        length_ray=zeros(len(vel))
        for klayer in range(len(vel)):
            if zshallow>depth1 and zdeep<depth2: #both points in same layer
                length_ray[klayer]=abs(zshallow-zdeep)/sin(theta) 
            elif zshallow>depth1 and zshallow<depth2: #This is the top
                length_ray[klayer]=abs(depth2-zshallow)/sin(theta)
            elif zdeep>depth1 and zdeep<depth2: #This is the bottom
                length_ray[klayer]=abs(depth1-zdeep)/sin(theta)
            elif depth1>zshallow and depth2<zdeep: #Use full layer thickness for ray path length
                length_ray[klayer]=delta_ray[klayer]
            else: #Some other layer, do nothing
                pass
            #Update reference depths
            if klayer<len(vel)-1: #last layer:
                depth1=depth2
                depth2=depth2+vel[klayer+1,0]
            else:
                depth1=depth2
                depth2=1e6
        
        
        #Now multiply ray path length times rupture velocity
        ray_times=length_ray/(vel[:,1]*rupture_multiplier)
        t_onset[kfault]=ray_times.sum()  
        length2fault[kfault]=(ray_times*vel[:,1]*rupture_multiplier).sum()
        
    #Now perturb onset times according to Graves-Pitarka eq 5 and 6 (assumes 1:1 corelation with slip)
    delta_t0=((M0*1e7)**(1./3))*1.8e-9
    
    #GP 2015 extra perturbation to destroy the 1:1 correlation with slip
    rand_numb=randn()
    delta_t=delta_t0*exp(sigma_rise_time*rand_numb)
    
    #Now apply total perturbation
    slip_average=slip.mean()
    i=where(slip>0.05*slip_average)[0] #perturbation is applied only to subfaults with significant slip
    perturbation=(log10(slip)-log10(slip_average))/(log10(slip.max())-log10(slip_average))
    t_onset_final=t_onset.copy()
    t_onset_final[i]=t_onset[i]-delta_t*perturbation[i]
    
    #Check for negative times
    i=where(t_onset_final<0)[0]
    t_onset_final[i]=t_onset[i]
    
    #print([M0,':',np.sum(np.sum(t_onset_final))])

    
    #Check for subfaults that would create P-waves arriving *earlier* than those from hypocenter, and redo those stochastic perturbations until they're all > hypo P-wave time
    for kfault in range(len(slip)):
        D,az,baz=gps2dist_azimuth(hypocenter[1],hypocenter[0],fault_array[kfault,2],fault_array[kfault,1])
        D=D/1000
        D_deg=kilometer2degrees(D)
        Ppaths=velmod.get_ray_paths(hypocenter[2],D_deg,phase_list=['P','p'])
        
        try:
            Ptime=Ppaths[0].time
        except:
            print('Ptime ERROR')
            print('   Target is:')
            print('     '+str(fault_array[kfault,0:4]))
            print('   Hypo is:')
            print('     '+str(hypocenter))
        attempt=1
        while (t_onset_final[kfault] < Ptime) & (attempt < 500):
            #if attempt == 1:
                #print("Trying at subfault " + str(kfault)
            #print(".....kfault " + str(kfault) + " is too early... recalculating!"
            rand_numb=randn()
            delta_t=delta_t0*exp(sigma_rise_time*rand_numb*(1./attempt))
            perturbation=(log10(slip[kfault])-log10(slip_average))/(log10(slip.max())-log10(slip_average))
            t_onset_final[kfault]=t_onset[kfault]-delta_t*perturbation
            attempt+=1
            if attempt == 500:
                t_onset_final[kfault]=t_onset[kfault]
                #print("...Failed at subfault " + str(kfault)
            #GP 2015 extra perturbation to destroy the 1:1 correlation with slip
            t_onset_final[kfault]=t_onset[kfault]

    #Ensure hypocenter onset time is zero
    t_onset_final[i_hypo]=0

    #Reassign subfaults within the "nucleation zone" to be their original, unperturbed onsets    
    #nu=some-relation-to-M0
#    i=where(t_onset<5.5)[0]
#    t_nucleation_edge=max(t_onset[i])
#    t_onset_final=t_onset_final+t_nucleation_edge
#    t_onset_final[i]=t_onset[i]
    #Check for nan times
    i=where(isnan(t_onset_final)==True)[0]
    t_onset_final[i]=0

    return t_onset_final,length2fault  
                

def get_centroid(fault):
    '''
    Get moment based centroid
    '''
    
    #Get slip at subfaults
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    #Get moment at subfaults
    moment=slip*fault[:,10]*fault[:,11]*fault[:,13]
    #centroid coordinates
    centroid_lon=sum(fault[:,1]*moment)/moment.sum()
    centroid_lat=sum(fault[:,2]*moment)/moment.sum()
    centroid_z=sum(fault[:,3]*moment)/moment.sum()

    return centroid_lon,centroid_lat,centroid_z             

                                                
def get_stochastic_rake(rake,Nsamples,sigma_rake=10,max_variation=45):
    '''
    Get stochastic rake
    '''
    
    from numpy.random import randn
    from numpy import where
    
    #Limits
    max_rake=rake+max_variation
    min_rake=rake-max_variation
    
    #Get stochastic rake
    stoc_rake=sigma_rake * randn(Nsamples) + rake
    
    #make sure we don't exceed the limits
    i=where(stoc_rake>max_rake)[0]
    stoc_rake[i]=max_rake
    i=where(stoc_rake<min_rake)[0]
    stoc_rake[i]=min_rake
    
    return stoc_rake                                                                                                
   
def write_all_event_summary(home,project_name,run_name):
    '''
    Write a sumamry file with Mw, max slip, rise_time,onset time for all events
    '''
    
    from glob import glob
    from string import replace
    from numpy import array,genfromtxt,sqrt,zeros,savetxt
    
    #How many events?
    ruptures=sorted(glob(home+project_name+'/output/ruptures/*.rupt'))
    logs=sorted(glob(home+project_name+'/output/ruptures/*.log'))
    
    #where does this go?
    fout=home+project_name+'/output/ruptures/_kin_summary.txt'
    out=zeros((len(logs),7))
    
    for k in range(len(logs)):
        
        print(k)
        
        #Get info about fault
        f=open(logs[k],'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Hypocenter (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])  
        f.close()        
        
        #Get peak quantities
        f=genfromtxt(ruptures[k])
        peak_slip=max(sqrt(f[:,8]**2+f[:,9]**2))
        peak_rise=max(f[:,7])
        peak_onset=max(f[:,12])
        
        out[k,0:3]=hypo
        out[k,3]=Mw
        out[k,4]=peak_slip
        out[k,5]=peak_rise
        out[k,6]=peak_onset
        
    savetxt(fout,out,fmt='%8.2f')
    
    
    
def build_TauPyModel(home,project_name,vel_mod_file,background_model='PREM'):
    '''
    This function will take the structure from the .mod file
    and paste it on top of a pre computed mantle structure such as PREM.
    This assumes that the .mod file provided by the user ends with a 0 thickness 
    layer on the MANTLE side of the Moho
    '''
    
    from numpy import genfromtxt
    from os import environ,path
    from obspy.taup import taup_create, TauPyModel
    
    #mudpy source folder
    
    #load user specified .mod infromation
    structure = genfromtxt(vel_mod_file)
    
    #load background velocity structure
    if background_model=='PREM':
        
        bg_model_file=environ['MUD']+'/src/aux/prem.nd'
        
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
        taup_create.build_taup_model(taup_in,output_folder=taup_out)
        
    else: #To be done later (ha)
        print('ERROR: That background velocity model does not exist')
        
    
    
    
    
def write_rupt_list(home,project_name,run_name,target_Mw,Nrealizations,ncpus):
    '''
    Writes ruptures.list file, note this does not check whether all ruptures
    were succesfully generated, it simply assumes that there were len(target_Mw)
    x Nrealizations ruptures made
    '''
   
    from numpy import ceil
    
    f=open(home+project_name+'/data/ruptures.list','w')
    
    #Check number of realizations based on ncpus
    if ncpus==1: #serial
        Nruptures = Nrealizations * len(target_Mw)
        
    else:
        #calculate number of realizations per CPU:
        Nrealizations_parallel=int(ceil(float(Nrealizations)/float(ncpus)))
        Nruptures = (Nrealizations_parallel*ncpus) * len(target_Mw)
    
    for krup in range(Nruptures):
    
        run_number = str(krup).rjust(6,'0')
        rupture = run_name+'.'+run_number+'.rupt\n'
        f.write(rupture)            

    f.close()

    
        
    

def generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
		load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,
		Lstrike,num_modes,Nrealizations,rake,rise_time,rise_time_depths,time_epi,
		max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,
		force_magnitude=False,force_area=False,mean_slip_name=None,hypocenter=None,
		slip_tol=1e-2,force_hypocenter=False,no_random=False,use_hypo_fraction=True,
		shear_wave_fraction_shallow=0.49,shear_wave_fraction_deep=0.8,max_slip_rule=True):
    '''
    Set up rupture generation-- use ncpus if available
    '''
    
    # shear_wave_fraction_shallow=1/60/60/24*3
    # shear_wave_fraction_deep=1/60/60/24*3
    
    #Things that need to be done before deciding whether to ship the work off to the 
    #parallel or serial functions
    
    #Need to make tauPy file
    vel_mod_file=home+project_name+'/structure/'+model_name
    #Get TauPyModel
    build_TauPyModel(home,project_name,vel_mod_file,background_model='PREM')

    #Write ruptures.list file
    write_rupt_list(home,project_name,run_name,target_Mw,Nrealizations,ncpus)
    
    
    #Generate rupture models
    run_generate_ruptures_parallel(home,project_name,run_name,fault_name,slab_name,mesh_name,
    load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,
    Lstrike,num_modes,Nrealizations,rake,rise_time,rise_time_depths,time_epi,
    max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,
    force_magnitude,force_area,mean_slip_name,hypocenter,slip_tol,force_hypocenter,
    no_random,use_hypo_fraction,shear_wave_fraction_shallow,shear_wave_fraction_deep,max_slip_rule)




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
def run_generate_ruptures_parallel(home,project_name,run_name,fault_name,slab_name,mesh_name,
        load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,
        Lstrike,num_modes,Nrealizations,rake,rise_time,rise_time_depths,time_epi,
        max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,
        force_magnitude,force_area,mean_slip_name,hypocenter,slip_tol,force_hypocenter,
        no_random,use_hypo_fraction,shear_wave_fraction_shallow,shear_wave_fraction_deep,max_slip_rule):
    
    from numpy import ceil
    from os import environ
    import subprocess
    from shlex import split
    
    #calculate number of realizations per CPU:
    Nrealizations_parallel=int(ceil(float(Nrealizations)/float(ncpus)))
    if (Nrealizations_parallel*ncpus > Nrealizations):
        print("Extra CPUS-- have " + str(Nrealizations_parallel*ncpus-Nrealizations) + " free ruptures!!")
    rise_time_depths0=rise_time_depths[0]
    rise_time_depths1=rise_time_depths[1]
    tMw=target_Mw[0]
    for r in range(len(target_Mw)-1):
        #tMw.append(target_Mw[r])
        tMw=str(tMw)+','+str(target_Mw[r+1])
    #Make mpi system call
    print("MPI: Starting " + str(Nrealizations_parallel*ncpus) + " FakeQuakes Rupture Generations on ", ncpus, "CPUs")
    mud_source=environ['MUD']+'/src/python/mudpy/'

    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'generate_ruptures_parallel.py run_parallel_generate_ruptures '+home+' '+project_name+' '+run_name+' '+fault_name+' '+str(slab_name)+' '+str(mesh_name)+' '+str(load_distances)+' '+distances_name+' '+UTM_zone+' '+str(tMw)+' '+model_name+' '+str(hurst)+' '+Ldip+' '+Lstrike+' '+str(num_modes)+' '+str(Nrealizations_parallel)+' '+str(rake)+' '+str(rise_time)+' '+str(rise_time_depths0)+' '+str(rise_time_depths1)+' '+str(time_epi)+' '+str(max_slip)+' '+source_time_function+' '+str(lognormal)+' '+str(slip_standard_deviation)+' '+scaling_law+' '+str(ncpus)+' '+str(force_magnitude)+' '+str(force_area)+' '+str(mean_slip_name)+' "'+str(hypocenter)+'" '+str(slip_tol)+' '+str(force_hypocenter)+' '+str(no_random)+' '+str(use_hypo_fraction)+' '+str(shear_wave_fraction_shallow)+' '+str(shear_wave_fraction_deep)+' '+str(max_slip_rule)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()

    
    
    
def run_generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
        load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,
        Lstrike,num_modes,Nrealizations,rake,rise_time,rise_time_depths,time_epi,
        max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,
        force_magnitude,force_area,mean_slip_name,hypocenter,slip_tol,force_hypocenter,
        no_random,shypo,use_hypo_fraction,shear_wave_fraction,max_slip_rule):
    
    '''
    Depending on user selected flags parse the work out to different functions
    '''
    
    from numpy import load,save,genfromtxt,log10,cos,sin,deg2rad,savetxt,zeros,where
    from time import gmtime, strftime
    from obspy.taup import TauPyModel


    #Should I calculate or load the distances?
    if load_distances==1:  
        Dstrike=load(home+project_name+'/data/distances/'+distances_name+'.strike.npy')
        Ddip=load(home+project_name+'/data/distances/'+distances_name+'.dip.npy')
    else:
        Dstrike,Ddip=subfault_distances_3D(home,project_name,fault_name,slab_name,UTM_zone)
        save(home+project_name+'/data/distances/'+distances_name+'.strike.npy',Dstrike)
        save(home+project_name+'/data/distances/'+distances_name+'.dip.npy',Ddip)
    

    #Read fault and prepare output variable
    whole_fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    
    #Get structure model
    vel_mod_file=home+project_name+'/structure/'+model_name
    
    #Get TauPyModel
    velmod = TauPyModel(model=home+project_name+'/structure/'+model_name.split('.')[0])


    #Now loop over the number of realizations
    realization=0
    print('Generating rupture scenarios')
    for kmag in range(len(target_Mw)):
        print('... Calculating ruptures for target magnitude Mw = '+str(target_Mw[kmag]))
        for kfault in range(Nrealizations):
            if kfault%1==0:
                print('... ... working on rupture '+str(kfault)+' of '+str(Nrealizations))
            
            #Prepare output
            fault_out=zeros((len(whole_fault),14))
            fault_out[:,0:8]=whole_fault[:,0:8]
            fault_out[:,10:12]=whole_fault[:,8:]   
            
            #Sucess criterion
            success=False
            while success==False:
                #Select only a subset of the faults based on magnitude scaling
                current_target_Mw=target_Mw[kmag]
                ifaults,hypo_fault,Lmax,Wmax,Leff,Weff=select_faults(whole_fault,Dstrike,Ddip,current_target_Mw,
                            num_modes,scaling_law,force_area,no_shallow_epi=False,
                            no_random=no_random,subfault_hypocenter=shypo,use_hypo_fraction=use_hypo_fraction)
                
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
                elif Ldip=='auto': #MH2019
                    Ld=6.8+0.4*Weff
                else:
                    Ld=Ldip
                
                #Get the mean uniform slip for the target magnitude
                if mean_slip_name==None:
                    mean_slip,mu=get_mean_slip(target_Mw[kmag],fault_array,vel_mod_file)
                else:
                    foo,mu=get_mean_slip(target_Mw[kmag],fault_array,vel_mod_file)
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
                C=vonKarman_correlation(Dstrike_selected,Ddip_selected,Ls,Ld,hurst)
                
                # Lognormal or not?
                if lognormal==False:
                    #Get covariance matrix
                    C_nonlog=get_covariance(mean_slip,C,target_Mw[kmag],fault_array,vel_mod_file,slip_standard_deviation) 
                    #Get eigen values and eigenvectors
                    eigenvals,V=get_eigen(C_nonlog)
                    #Generate fake slip pattern
                    rejected=True
                    while rejected==True:
#                        slip_unrectified,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=kfault)
                        slip_unrectified,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=None)
                        slip,rejected,percent_negative=rectify_slip(slip_unrectified,percent_reject=13)
                        if rejected==True:
                            print('... ... ... negative slip threshold exceeeded with %d%% negative slip. Recomputing...' % (percent_negative))
                else:
                    #Get lognormal values
                    C_log,mean_slip_log=get_lognormal(mean_slip,C,target_Mw[kmag],fault_array,vel_mod_file,slip_standard_deviation)               
                    #Get eigen values and eigenvectors
                    eigenvals,V=get_eigen(C_log)
                    #Generate fake slip pattern
#                    slip,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=kfault)
                    slip,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=None)
            
                #Slip pattern sucessfully made, moving on.
                #Rigidities
                foo,mu=get_mean_slip(target_Mw[kmag],whole_fault,vel_mod_file)
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
            stoc_rake=get_stochastic_rake(rake,len(slip))
            
            #Place slip values in output variable
            fault_out[ifaults,8]=slip*cos(deg2rad(stoc_rake))
            fault_out[ifaults,9]=slip*sin(deg2rad(stoc_rake))
            
#            #Move hypocenter to somewhere with a susbtantial fraction of peak slip
#            slip_fraction=0.25
#            islip=where(slip>slip.max()*slip_fraction)[0]
#            shuffle(islip) #randomize
#            hypo_fault=ifaults[islip[0]] #select first from randomized vector
            
            #Calculate and scale rise times
            rise_times=get_rise_times(M0,slip,fault_array,rise_time_depths,stoc_rake,rise_time)
            
            #Place rise_times in output variable
            fault_out[:,7]=0
            fault_out[ifaults,7]=rise_times
            
            #Calculate rupture onset times
            if force_hypocenter==False: #Use random hypo, otehrwise force hypo to user specified
                hypocenter=whole_fault[hypo_fault,1:4]
            
            t_onset=get_rupture_onset(home,project_name,slip,fault_array,model_name,hypocenter,
                                      rise_time_depths,M0,velmod,shear_wave_fraction)
            fault_out[:,12]=0
            fault_out[ifaults,12]=t_onset
            
            #Calculate location of moment centroid
            centroid_lon,centroid_lat,centroid_z=get_centroid(fault_out)
            
            #Write to file
            run_number=str(realization).rjust(6,'0')
            outfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.rupt'
            savetxt(outfile,fault_out,fmt='%d\t%10.6f\t%10.6f\t%8.4f\t%7.2f\t%7.2f\t%4.1f\t%5.2f\t%5.2f\t%5.2f\t%10.2f\t%10.2f\t%5.2f\t%.6e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
            
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
            f.write('Source time function type: %s' % source_time_function)
            f.close()
                        
            realization+=1

    
