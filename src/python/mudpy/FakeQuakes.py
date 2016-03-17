'''
D.Melgar 03/2016

Generate syntehtic slip distributions using the K-L expansion method.
This is modified from RJ Leveque's KL2d_vonKarman notebook to work in 
a similar way to the MudPy ivnersions
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
        clob=raw_input('Project directory exists, clobber (y/n)?')
        if clob is'y' or clob is 'Y': #Clobber baby
            clob=raw_input('This will delete everything in this project directory, so, take a minute, think about it: clobber (y/n)?')
            if clob is 'y' or clob is 'Y':
                rmtree(proj_dir)
            else: #Leave direcory alone
                print 'Phew, almost shot yourself in the foot there didn\'t you?'
        else: #Leave direcory alone
            print 'Phew, almost shot yourself in the foot there didn\'t you?'
    if clob is 'y' or clob is 'Y':
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




def subfault_distances_3D(home,project_name,fault_name,slab_name,projection_zone):
    """
    Estimate the distance between subfaults i and j for every pair in the list
    fault.subfaults. For a #D fault geometry

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

    from numpy import sqrt,sin,cos,deg2rad,zeros,meshgrid,linspace,where,c_,unravel_index,sort,diff,genfromtxt,sign
    from matplotlib.mlab import griddata
    from matplotlib import pyplot as plt
    from gmsh_tools import llz2utm
    from scipy.spatial.distance import cdist
    
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
    i=where(down_dip>360)[0]
    down_dip[i]=down_dip[i]-360
    
    #Convert slab1.0 to local UTM coordinates
    slab_x,slab_y=llz2utm(slab_model[:,0],slab_model[:,1],projection_zone)
    slab_x,slab_y = slab_x/1000,slab_y/1000
    slab_z=-slab_model[:,2]
    
    #Convert faul centroid coordinates to local UTM
    fault_x,fault_y=llz2utm(fault[:,1],fault[:,2],projection_zone)
    fault_x,fault_y = fault_x/1000,fault_y/1000
    fault_z=fault[:,3]    
    
    # grid Slab1.0 for making depth contours to be then used for along-strike distance calculation
    ngrid_pts=500
    X=linspace(slab_x.min(),slab_x.max(),ngrid_pts)
    Y=linspace(slab_y.min(),slab_y.max(),ngrid_pts)
    X,Y = meshgrid(X,Y)
    Z = griddata(slab_x, slab_y, slab_z, X, Y,interp='linear')
    
    # X, Y and Z are matrices with the grid info, now create one contour at each subfault centroid depth
    contour_levels=fault[:,3]
    all_contours=plt.contour(X,Y,Z,levels=contour_levels)
    
    # x-coordinates for down_dip line
    x_range=slab_x.max()-slab_x.min()
    x_down_dip=linspace(-x_range/2,x_range/2,200)
    
    #Loop over number of subfaults, we want the distance from i-th fault to all other (j) subfaults
    print 'Getting inter-fault distances'
    for i in range(len(fault)):
        if i%10==0:
            print '... working on subfault '+str(i)+' of '+str(len(fault))
        #Current fault
        xi = fault_x[i]
        yi = fault_y[i]
        zi = fault_z[i]
        
        #Get contour at depth of current subfault
        contour=all_contours.collections[i].get_paths()[0].vertices
        
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
            y_down_dip=x_down_dip*cos(deg2rad(down_dip[j]))
            
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
    
    #Print to keep track of things
    print "Correlation lengths: Lstrike = %g, Ldip = %g" % (Lstrike,Ldip)
    print "Hurst exponent = %4.2f" % hurst
    r = sqrt((Dstrike/Lstrike)**2 + (Ddip/Ldip)**2)
    C = vonKarman(r,hurst)
    
    return C


def get_lognormal(C,target_Mw,fault_array,vel_mod,alpha=0.8):
    '''
    Exponentiate to get lognormal correlation
    '''
    
    from numpy import log,diag
    
    #First get correct mean slip
    mean_slip,mu=get_mean_slip(target_Mw,fault_array,vel_mod)
    #Get sigma as a fraction of mean slip
    sigma_slip=alpha*mean_slip
    #Now generate desired correlation
    Cov = sigma_slip * (C*sigma_slip).T
    Cov_g = log((sigma_slip/mean_slip) * (C*(sigma_slip/mean_slip)).T + 1.)
    mean_slip_g = log(mean_slip) - diag(Cov_g)/2.
    
    return Cov_g,mean_slip_g
       
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
    
    
def make_KL_slip(fault,num_modes,eigenvals,V,mean_slip):
    '''
    Make slip map using num_modes
    '''
    from numpy import sqrt,exp
    from numpy.random import randn
    
    #Generate random numbers
    z = randn(num_modes) 
    KL_slip = mean_slip.copy()  # start with the mean slip
    # add in the terms in the K-L expansion:
    for k in range(len(z)):
        KL_slip += z[k] * sqrt(eigenvals[k]) * V[:,k]
    # exponentiate for lognormal:
    KL_slip = exp(KL_slip)
    #close("all");scatter(fault_x,fault_y,c=KL_slip,lw=0);axis('equal');colorbar()

    return KL_slip
    

def plot_KLslip(home,project_name,run_name,run_number,fault,mesh_name,Mw,hypo_fault,fudge=0.3):
    '''
    Make simple plot of slip model
    '''
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib import pyplot as plt
    from matplotlib import cm,colors
    from numpy import r_,genfromtxt,reshape
    from numpy.random import randn
    
    #Read mesh
    mesh=genfromtxt(home+project_name+'/data/model_info/'+mesh_name)
    
    #Colormap
    cdict = {'red': ((0., 1, 1),
                 (0.10, 1, 1),
                 (0.20, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.10, 1, 1),
                   (0.20, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.15, 1, 1),
                  (0.20, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
    whitejet = colors.LinearSegmentedColormap('whitejet',cdict,256)
    
    #Plot
    fig, ax = plt.subplots(figsize=(3,8))
    ax.set_xlim([fault[:,1].min()-fudge,fault[:,1].max()+fudge])
    ax.set_ylim([fault[:,2].min()-fudge,fault[:,2].max()+fudge])

    #Make patches
    patches = []
    for k in range(len(fault)):
        coordinates=reshape(r_[mesh[k,4:6],mesh[k,7:9],mesh[k,10:12],mesh[k,4:6]],(4,2))
        subfault = Polygon(coordinates, True)
        patches.append(subfault)
        
    p = PatchCollection(patches,cmap=whitejet,lw=0.2)
    colors = fault[:,9]
    p.set_array(colors)
    plt.xticks(rotation=70)
    ax.add_collection(p)
    plt.scatter(fault[hypo_fault,1],fault[hypo_fault,2],marker='s',lw=0.5,c='m',s=60)
    plt.title('Mw = %.2f' % (Mw))
    cbar_ax = fig.add_axes([0.73, 0.2, 0.05, 0.6])
    plt.colorbar(p, cax=cbar_ax, label="Slip (m)")
    plt.subplots_adjust(right=0.75)
    plt.savefig(home+project_name+'/plots/'+run_name+'.'+run_number+'.png')
    plt.close()


def select_faults(whole_fault,Dstrike,Ddip,target_Mw,buffer_factor):
    '''
    Select a random fault to be the hypocenter then based on scaling laws and a 
    target magnitude select only faults within the expected area plus some 
    buffer factor
    '''
    
    from numpy.random import randint
    from numpy import array,where
    
    done=False
    while not done:
        #Determine length and width from scaling laws
        length=10**(-2.37+0.57*target_Mw)
        length=length+length*buffer_factor
        width=10**(-1.86+0.46*target_Mw)
        width=width+width*buffer_factor
        
        #Select random subfault as hypocenter
        hypo_fault=randint(0,len(whole_fault)-1)

        #Get max/min distances from hypocenter to all faults
        dstrike_max=Dstrike[:,hypo_fault].max()
        dstrike_min=Dstrike[:,hypo_fault].min()
        dstrike_range=dstrike_max-dstrike_min
        ddip_max=Ddip[:,hypo_fault].max()
        ddip_min=Ddip[:,hypo_fault].min()
        ddip_range=ddip_max-ddip_min
        
        #Work on strike first
        strike_bounds=array([-length/2,length/2])
        if strike_bounds[0]<dstrike_min:#Length is outside domain
            strike_bounds[1]=strike_bounds[1]+abs(dstrike_min-strike_bounds[0])
            strike_bounds[0]=dstrike_min
        if strike_bounds[1]>dstrike_max:#Length is outside domain
            strike_bounds[0]=strike_bounds[0]-abs(dstrike_max-strike_bounds[1])
            strike_bounds[1]=dstrike_max
            
        #Now get dip ranges
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
        #From within the selected faults randomly select the hypocenter
        if len(selected_faults)<2:
            done=False
        else:
            i=randint(0,len(selected_faults)-1)
            hypo_fault=selected_faults[i]
            break
    return selected_faults,hypo_fault
    
def generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
    load_distances,distances_name,UTM_zone,target_Mw,model_name,hurst,Ldip,Lstrike,
    num_modes,Nrealizations,rake,buffer_factor):
    
    '''
    Depending on user selected flags parse the work out to different functions
    '''
    
    from numpy import load,save,genfromtxt,log10,cos,sin,deg2rad,savetxt,zeros,sqrt
    from string import rjust
    
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
    vel_mod=home+project_name+'/structure/'+model_name
    
    #Determine correlation lengths
    if Lstrike=='auto': #Use scaling
        Ls=10**(-2.43+0.49*target_Mw)
    else:
        Ls=Lstrike
    if Ldip=='auto': #Use scaling
        Ld=10**(-1.79+0.38*target_Mw)
    else:
        Ld=Ldip

    #Now loop over the number of realizations
    print 'Generating rupture scenarios'
    for kfault in range(Nrealizations):
        if kfault%10==0:
            print '... working on rupture '+str(kfault)+' of '+str(Nrealizations)
        
        #Prepare output
        fault_out=zeros((len(whole_fault),14))
        fault_out[:,0:8]=whole_fault[:,0:8]
        fault_out[:,10:12]=whole_fault[:,8:]   
        
        #Select only a subset of the faults based on magnitude scaling
        ifaults,hypo_fault=select_faults(whole_fault,Dstrike,Ddip,target_Mw,buffer_factor)
        fault_array=whole_fault[ifaults,:]
        Dstrike_selected=Dstrike[ifaults,:][:,ifaults]
        Ddip_selected=Ddip[ifaults,:][:,ifaults]
        
        #Get the mean uniform slip for the target magnitude
        mean_slip,mu=get_mean_slip(target_Mw,fault_array,vel_mod)
        
        #Get correlation matrix
        C=vonKarman_correlation(Dstrike_selected,Ddip_selected,Ls,Ld,hurst)
        
        #Get lognormal valcues
        C_log,mean_slip_log=get_lognormal(C,target_Mw,fault_array,vel_mod)
        
        #Get eigen values and eigenvectors
        eigenvals,V=get_eigen(C_log)
        
        #Generate fake slip pattern
        slip=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log)
        
        #Place in output variable
        fault_out[ifaults,8]=slip*cos(deg2rad(rake))
        fault_out[ifaults,9]=slip*sin(deg2rad(rake))
        
        #Rigidities
        foo,mu=get_mean_slip(target_Mw,whole_fault,vel_mod)
        fault_out[:,13]=mu
        
        #Calculate moment and magnitude of fake slip pattern
        M0=sum(sqrt(fault_out[:,8]**2+fault_out[:,9]**2)*fault_out[:,10]*fault_out[:,11]*mu)
        Mw=(2./3)*(log10(M0)-9.1)
        
        #Write to file
        run_number=rjust(str(kfault),4,'0')
        outfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.rupt'
        savetxt(outfile,fault_out,fmt='%d\t%10.6f\t%10.6f\t%7.2f\t%7.2f\t%7.2f\t%4.1f\t%5.2f\t%5.2f\t%5.2f\t%10.2f\t%10.2f\t%5.2f\t%.6e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
        
        #Make plots
        plot_KLslip(home,project_name,run_name,run_number,fault_out,mesh_name,Mw,hypo_fault,fudge=0.3)
        
        #Write log file
        logfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(logfile,'w')
        f.write('Project name: '+project_name+'\n')
        f.write('Run name: '+run_name+'\n')
        f.write('Run number: '+run_number+'\n')
        f.write('Velocity model: '+model_name+'\n')
        f.write('No. of KL modes: '+str(num_modes)+'\n')
        f.write('Hurst exponent: '+str(hurst)+'\n')
        f.write('Lstrike: '+str(Ls)+'km\n')
        f.write('Ldip: '+str(Ld)+'km\n')
        f.write('Target magnitude: Mw '+str(target_Mw)+'\n')
        f.write('Actual magnitude: Mw '+str(Mw))
        f.close()
        