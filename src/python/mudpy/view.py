'''
D.Melgar
04/2014

Some routines to make quick and not so quick plots of the forward modeling and 
inversion results
'''

import matplotlib
from matplotlib import cm
from mudpy.gmttools import gmtColormap

#Create default colormap for slip inversions
#cdict = {'red': ((0., 1, 1),
#                 (0.05, 1, 1),
#                 (0.11, 0, 0),
#                 (0.66, 1, 1),
#                 (0.89, 1, 1),
#                 (1, 0.5, 0.5)),
#         'green': ((0., 1, 1),
#                   (0.05, 1, 1),
#                   (0.11, 0, 0),
#                   (0.375, 1, 1),
#                   (0.64, 1, 1),
#                   (0.91, 0, 0),
#                   (1, 0, 0)),
#         'blue': ((0., 1, 1),
#                  (0.05 1, 1),
#                  (0.11, 1, 1),
#                  (0.34, 1, 1),
#                  (0.65, 0, 0),
#                  (1, 0, 0))}
cdict = {'red': ((0., 1, 1),
                 (0.03, 1, 1),
                 (0.20, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.03, 1, 1),
                   (0.20, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.08, 1, 1),
                  (0.20, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

#cdict = {'red': ((0., 1, 1),
#                 (0.03, 1, 1),
#                 (0.20, 0, 0),
#                 (0.66, 1, 1),
#                 (0.89, 1, 1),
#                 (1, 0.5, 0.5)),
#         'green': ((0., 1, 1),
#                   (0.03, 1, 1),
#                   (0.20, 0, 0),
#                   (0.375, 1, 1),
#                   (0.64, 1, 1),
#                   (0.91, 0, 0),
#                   (1, 0, 0)),
#         'blue': ((0., 1, 1),
#                  (0.045, 1, 1),
#                  (0.20, 1, 1),
#                  (0.34, 1, 1),
#                  (0.65, 0, 0),
#                  (1, 0, 0))}
#whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)


#whitejet=gmtColormap(u'/Users/dmelgarm/code/python/cpt/color_linear.cpt')

pqlx_dict={'blue': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  1.  ,  1.  ),
            ( 0.35 ,  1.  ,  1.  ),
            ( 0.6 ,  1.  ,  1.  ),
            ( 0.7 ,  0.  ,  0.  ),
            ( 0.89 ,  0.  ,  0.  ),
            ( 1.  ,  0.  ,  0.  )), 
            'green': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  0.  ,  0.  ),
            ( 0.35 ,  0.  ,  0. ),
            ( 0.6 ,  1.  ,  1.  ),
            ( 0.7 ,  1.  ,  1.  ),
            ( 0.89 ,  1.  ,  1.  ),
            ( 1.  ,  0.  ,  0. )), 
            'red': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  1.  ,  1.  ),
            ( 0.35 ,  0.  ,  0.  ),
            ( 0.6 ,  0.  ,  0.  ),
            ( 0.7 ,  0.  ,  0.  ),
            ( 0.89 ,  1.  ,  1.  ),
            ( 1.  ,  1.  ,  1.  ))}
            

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)


def quick_model(rupt,s=5,slip_percent=0,facecolor='r'):
    '''
    Quick and dirty plot of a .rupt file. Shows map view of slip
    
    Parameters:
            rupt: string
            The absolute path to a .inv or .rupt file
            
    Example:
        view.quick_model('/foo/bar/output/inverse_models/models/inv_result.inv')
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import matplotlib.pyplot as plt
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get other parameters
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    strike=f[0:len(unum),4]
    #Get projection of rake vector
    x,y=slip2geo(ss,ds,strike)
    
    #keep only appropriate slip
    i=where(slip > slip_percent*slip.max())
    slip=slip[i]
    lon=lon[i]
    lat=lat[i]
    x=x[i]
    y=y[i]
    
    #Plot
    plt.figure(figsize=(5.2,10))
    plt.scatter(lon,lat,marker='o',c=slip,s=s,cmap=whitejet,vmin=0,facecolor=facecolor)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.quiver(lon,lat,x,y,color='green',width=0.0013)
    plt.grid()
    #plt.title(rupt)
    #plt.title(rupt[58:85])
    #plt.savefig(rupt[0:49]+'figures/'+rupt[58:85]+'.png')
    #plt.close()
    plt.show()


def plot_onsettime(rupt,s=5):
    '''
    Quick and dirty plot of a .rupt file. Shows map view of slip
    
    Parameters:
            rupt: string
            The absolute path to a .inv or .rupt file
            
    Example:
        view.quick_model('/foo/bar/output/inverse_models/models/inv_result.inv')
    '''
    
    from numpy import genfromtxt,unique,where,zeros,array
    import matplotlib.pyplot as plt
    from string import replace
    
    f=genfromtxt(rupt)
    onset=f[:,12]
    num=f[:,0]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    #Get other parameters
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    #Find hypocenter
    log_file=rupt[:-5]+'.log'
    j=open(log_file,'r')
    loop_go=True
    while loop_go:
        line=j.readline()  
        if 'Hypocenter (lon,lat,z[km])' in line:
            q=replace(line.split(':')[-1],'(','')
            q=replace(q,')','')
            hypocenter=array(q.split(',')).astype('float')
            loop_go=False       
    #Plot
    plt.figure(figsize=(5.2,10))
    plt.scatter(hypocenter[0],hypocenter[1],marker='o',edgecolors='k',lw=2,s=5*s)
    plt.scatter(lon,lat,marker='o',c=onset,s=s,cmap=whitejet,vmax=5)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.xlim([hypocenter[0]-1,hypocenter[0]+1])
    plt.ylim([hypocenter[1]-1,hypocenter[1]+1])
    cb=plt.colorbar()
    cb.set_label('Onset Time (s)')
    plt.grid()
    plt.title(rupt[-26:-5])
    #plt.title(rupt[58:85])
    plt.savefig(rupt[:-35]+'figures/'+rupt[-26:-5]+'/'+rupt[-26:-5]+'_onset.png')
    plt.close()
    #plt.show()


    
def quick_static(gflist,datapath,scale=1):
    '''
    Make quick quiver plot of static field data
    
    Parameters:
        gflist: string
        Absolute path to GF_list file
        
        datapath: string
        Absolute path to data files
        
        scale: float, optional
        Scale value for quiver
        
        run_name: string, optional
        Run name of inversion
        
        run_number: string, optional
        Run number of inversion
        
    Examples:
        
        If plotting input static data
        >>> gflist=u'/Users/dmelgar/Slip_inv/Nepal/data/station_info/GPS.gflist'
        >>> datapath='/Users/dmelgar/Slip_inv/Nepal/data/statics/'
        >>> view.quick_static(gflist,datapath,scale=10)
        
    Notes:
        The value of scale works counterintuitvely. A larger value makes the 
        arrow lengths smaller and viceversa.
        
       This code only works if files are named according to the format sta.neu       
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,where,zeros,meshgrid,linspace

    
    GF=genfromtxt(gflist,usecols=3)
    sta=genfromtxt(gflist,usecols=0,dtype='S')
    lon=genfromtxt(gflist,usecols=1,dtype='f')
    lat=genfromtxt(gflist,usecols=2,dtype='f')
    #Read coseismcis
    i=where(GF!=0)[0]
    lon=lon[i]
    lat=lat[i]
    n=zeros(len(i))
    e=zeros(len(i))
    u=zeros(len(i))
    #Get data
    plt.figure()
    for k in range(len(i)):
        plt.annotate(sta[k],xy=(lon[k],lat[k]))
        neu=genfromtxt(datapath+sta[i[k]]+'.neu')
        n[k]=neu[0]
        e[k]=neu[1]
        u[k]=neu[2]          
    #Plot
    plt.quiver(lon,lat,e,n,scale=scale)
    plt.scatter(lon,lat,color='b')
    plt.grid()
    plt.title(datapath)
    plt.show()


def slip3D(rupt,marker_size=60,clims=None,plot_onset=False,cmap=whitejet,show=False):
    '''
    For complex fault geometries make a quick 3D plot of the rupture model
    
    Parameters:
            rupt: string
            The absolute path to a .inv or .rupt file
            
            marker_size: int, optional
            The size of the subfault markers, defaults to 60
            
    Example:
        >>> rupt='/Users/dmelgar/Slip_inv/Nepal/output/inverse_models/models/final.0000.inv'
        >>> view.slip3D(rupt,80)
    '''
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,zeros,unique,where

    #Parse rupture or inverse file
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get other parameters
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    depth=-f[0:len(unum),3]

    #get onsets
    onsets=f[0:len(unum),12]

    if plot_onset==False:
        plot_variable=slip
    else:
        plot_variable=onsets

    #Plot it
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(111, projection='3d')
    if clims==None:
        p=ax.scatter(lon, lat, depth, c=plot_variable,cmap=cmap, marker='o',s=marker_size,lw=0)
    else:
        p=ax.scatter(lon, lat, depth, c=plot_variable,cmap=cmap, marker='o',s=marker_size,vmin=clims[0],vmax=clims[1],lw=0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth (km)')
    cb=fig.colorbar(p)
    
    if plot_onset==False:
        cb.set_label('Slip (m)')
    else:
        cb.set_label('Onset time (s)')
    plt.subplots_adjust(left=0.1, bottom=0.1, right=1.0, top=0.9, wspace=0, hspace=0)
    plt.title(rupt)
    if show:
        plt.show()



def fence_slip(home,project_name,run_name,run_number,maxslip=None,UTM_zone='10S',elev=20,azimuth=None,fudge=10,fake_hypo=[0.1,0.1],
        borderwidth=0.5,figsize=(21,4),xtick=10,ytick=10,ztick=5,hypocenter=None):
    '''
    Make fence diagram of rupture file
    '''
    
    from numpy import genfromtxt,array,zeros,where
    from matplotlib import pyplot as plt
    from matplotlib import colors
    from pyproj import Proj
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from string import replace
    from matplotlib import ticker
    from matplotlib.ticker import MultipleLocator
    from mudpy.viewFQ import get_subfault_corners,corners2utm
    from mudpy import gmttools

    #Get rupture data
    fault_name=home+project_name+'/output/inverse_models/models/%s.%s.inv' % (run_name,run_number)
    gmttools.make_total_model(fault_name,thresh=0)
    fault=genfromtxt(home+project_name+'/output/inverse_models/models/%s.%s.inv.total' % (run_name,run_number))
    #Parse log file for hypocenter
    log_file=home+project_name+'/output/inverse_models/models/%s.%s.log' % (run_name,run_number)
    f=open(log_file,'r')
    loop_go=True
    while loop_go:
        line=f.readline()  
        if 'Mw' in line:
            Mw=float(line.split(':')[-1].split(' ')[-1])   
            loop_go=False
    f.close() 
    
    #get subfault corners
    corners=get_subfault_corners(fault)
    
    #Convert ot UTM (in km)
    corners=corners2utm(corners,UTM_zone=UTM_zone)
    
    #Get UTM coords of hypocenter
    P=Proj(proj='utm',ellps='WGS84',zone=UTM_zone)
    hypocenter[0],hypocenter[1]=P(hypocenter[0],hypocenter[1])
    hypocenter[0]/=1000
    hypocenter[1]/=1000
    
    #Make hypocenter the origin
    corners[:,0]-=hypocenter[0]
    corners[:,3]-=hypocenter[0]
    corners[:,6]-=hypocenter[0]
    corners[:,9]-=hypocenter[0]
    corners[:,1]-=hypocenter[1]
    corners[:,4]-=hypocenter[1]
    corners[:,7]-=hypocenter[1]
    corners[:,10]-=hypocenter[1]
    
    #Get mean strike for inital viewing angle
    strike=fault[:,4].mean()
    
    #Normalized slip
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    
    #Saturate to maxslip
    if maxslip!=None:
        imax=where(slip>maxslip)[0]
        slip[imax]=maxslip
    #normalize
    norm_slip=slip/slip.max()
    
    #Get colormaps
    pqlx = colors.LinearSegmentedColormap('pqlx',pqlx_dict,256)
    
    #Azimuth viewing angle
    if azimuth==None:
        azimuth=strike+90
    
    #Plot init, axes positions etc
    fig=plt.figure(figsize=figsize)
    
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_xlim([corners[:,0].min()-fudge,corners[:,0].max()+fudge])
    ax1.set_ylim([corners[:,1].min()-fudge,corners[:,1].max()+fudge])
    ax1.set_zlim([corners[:,2].min()-fudge/4,corners[:,2].max()+fudge/4])
    #Fenagle the axis ticks
    xmajorLocator = MultipleLocator(xtick)
    ymajorLocator = MultipleLocator(ytick)
    zmajorLocator = MultipleLocator(ztick)
    ax1.xaxis.set_major_locator(xmajorLocator)
    ax1.yaxis.set_major_locator(ymajorLocator)
    ax1.zaxis.set_major_locator(zmajorLocator)
    ax1.invert_zaxis()
    ax1.view_init(elev=elev, azim=azimuth)

    #Make one patch per subfault
    for ksub in range(len(corners)):
        vertices=[[tuple(corners[ksub,0:3]),tuple(corners[ksub,3:6]),tuple(corners[ksub,6:9]),tuple(corners[ksub,9:12])]]
        subfault=Poly3DCollection(vertices, linewidths=borderwidth)
        #subfault.set_color(pqlx(norm_slip[ksub]))
        subfault.set_color(pqlx(norm_slip[ksub]))
        subfault.set_linewidth(borderwidth)
        subfault.set_edgecolor('#505050')
        ax1.add_collection3d(subfault)
    
    #Hypocenter
    ax1.scatter(fake_hypo[0],fake_hypo[1],hypocenter[2],s=220,marker='*',c='#7CFC00')
    
    #Dummy mapable for colorbar
    s=plt.scatter(zeros(len(fault)),zeros(len(fault)),c=slip,cmap=pqlx,s=0.00001,lw=0)
    
    #Mke colorbar
    cb=plt.colorbar(s,shrink=0.9,pad=-0.07)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    cb.set_label('Slip (m)')
    
    #Labels n' stuff
    ax1.set_xlabel('\n\nEast (km)')
    ax1.set_ylabel('\n\nNorth (km)')
    ax1.set_zlabel('Depth (km)',rotation=90)
    #plt.title(home+project_name+'/output/ruptures/%s.%s.rupt' % (run_name,run_number))
    plt.title(run_name+' '+run_number+' Mw '+str(Mw))
    
    plt.show()





def plot_insar(home,project_name,GF_list,los_min,los_max):
    '''
    Plot the InSAR LOS data
    '''
    from numpy import genfromtxt,where,zeros
    from matplotlib import pyplot as plt
    from matplotlib import cm
    path=home+project_name+'/data/statics/'
    stations=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    i=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=7)
    i=where(i==1)[0]
    lon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=1)
    lat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=2)
    stations=stations[i]
    lon=lon[i]
    lat=lat[i]
    los=zeros((len(lon),4))
    for k in range(len(lon)):
        los[k,:]=genfromtxt(path+stations[k]+'.los')
    plt.figure()
    plt.subplot(221)
    plt.scatter(lon,lat,c=los[:,0],cmap=cm.jet,vmin=los_min,vmax=los_max,lw=0)
    plt.colorbar()
    plt.title('LOS(m)')
    f=genfromtxt('/Users/dmelgar/Slip_inv/Lefkada70/data/model_info/lefkada65.fault')
    #plt.scatter(f[:,1],f[:,2],marker='x')
    plt.subplot(222)
    plt.scatter(lon,lat,c=los[:,1],cmap=cm.jet,lw=0)
    plt.colorbar()
    plt.title('North') 
    plt.subplot(223)
    plt.scatter(lon,lat,c=los[:,2],cmap=cm.jet,lw=0)
    plt.colorbar()
    plt.title('East') 
    plt.subplot(224)
    plt.scatter(lon,lat,c=los[:,3],cmap=cm.jet,lw=0)
    plt.colorbar()
    plt.title('Up') 
    plt.show()


def tile_slip(rupt,nstrike,ndip,slip_bounds,geographic=False,epicenter=0,epicenter_line=0,
              thresh=0,xlims=[-100,100],ylims=[-100,100],fig_size=(10, 3),cmap=whitejet,
              afters=False,afters_file=None,histograms=False,size=250,vertices=None):
    '''
    Detailed plot of a forward model or inversion result file
    
    Parameters:
        rupt: string
        Absolute path to forward model (.rupt) or inversion result (.inv) file
        
        nstrike: int
        Number of along-strike subfaults in fault model
        
        ndip: int
        Number of down -dip subfaults in fault model
        
        slip_bounds: tuple, (slip_min(float),slip_max(float))
        Limits of slip used for colorscale
        
        geographic: boolean, optional
        Used to decide between plotting in geographic coordinates or subfault index number
        
        epicenter: numpy array, optional
        Epicentral coordinates as [longitude,latitude,depth(km)]
        
        epicenter_line: int, optional
        Only used if geographic=True, down-dip line number (starting from 1) where
        hypocenter is located
        
        thresh, float, optional
        Do not plot slip values below this value
        
    Example:
        >>> view.tile_slip('/Users/dmelgar/Slip_inv/Nepal/output/inverse_models/models/GPS.000.inv',
                20,15,(0,6.5),geographic=True,epicenter=epicenter,epicenter_line=8)
                
    Notes:
        If epicenter_line is not properly set and geographic=True then the aspect ratio
        of the subfaults will look wrong.
        
        This routine will plot the rake vectors determined from the input file.
        
        If saved as .pdf this makes a publication quality plot.
    '''
    
    from numpy import genfromtxt,unique,where,zeros,tile,sqrt
    import matplotlib.pyplot as plt
    from obspy.geodetics import gps2dist_azimuth
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Apply threshold
    ithresh=where(slip<thresh)[0]
    slip[ithresh]=0
    #Get unit rake vector
    rakess=ss/slip
    rakeds=ds/slip
    slip_min=slip_bounds[0]
    slip_max=slip_bounds[1]
    #Aftershocks
    if afters==True:
        af=genfromtxt(afters_file)
#        af=genfromtxt('/Users/dmelgarm/Turkey2020/catalogs/afters_10km.txt')
    #    af=genfromtxt('/Users/dmelgarm/Turkey2020/catalogs/pre_10km.txt')
        lon_afters=af[:,0]
        lat_afters=af[:,1]
        depth_afters=af[:,2]

    if geographic==True: #Get geographic coordinates to compute along strike and along dip distance
        lon=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,1] #Only compute line at the epicenter depth
        lat=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,2]    
        depth=-f[:,3]
        depth=depth[0:len(unum)]
        along_strike=zeros(nstrike)
        for k in range(len(lat)):
            out=gps2dist_azimuth(epicenter[1],epicenter[0],lat[k],lon[k])
            if lat[k]>epicenter[1]: #It's to the south
                #along_strike[k]=-out[0]/1000
                along_strike[k]=out[0]/1000
            else:
                #along_strike[k]=out[0]/1000
                along_strike[k]=-out[0]/1000
        #Now tile
        along_strike=tile(along_strike,ndip)
        
        
        #Process the aftershocks
        if afters==True:
            along_strike_afters=zeros(len(lon_afters))
            for k in range(len(lat_afters)):
                out=gps2dist_azimuth(epicenter[1],epicenter[0],lat_afters[k],lon_afters[k])
                if lat_afters[k]<epicenter[1]: #It's to the south
                    #along_strike_afters[k]=-out[0]/1000
                    along_strike_afters[k]=out[0]/1000
                else:
                    along_strike_afters[k]=-out[0]/1000
    
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1          
    #Plot
    if geographic==False:
        plt.figure()
        plt.scatter(istrike,idip,marker='o',c=slip,s=size,cmap=cmap,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
        plt.ylabel('Along-dip index')
        plt.xlabel('Along-strike index')
        plt.xlim(istrike.min()-1,istrike.max()+1)
        plt.ylim(idip.min()-1,idip.max()+1)
        #plt.quiver(istrike,idip,rakess,rakeds,color='green',width=0.005,scale=60.0)
        plt.axis('equal')
        plt.grid()    
        plt.title(rupt)
    else:
        rakess=rakess*slip
        rakeds=rakeds*slip
        plt.figure(num=None, figsize=fig_size, dpi=80)
        if vertices != None:
            plt.scatter(along_strike,depth,marker=vertices,linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=size,cmap=cmap,vmin=slip_min,vmax=slip_max)
        else:
            plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=size,cmap=cmap,vmin=slip_min,vmax=slip_max)
       
            #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=plt.cm.afmhot_r,vmin=slip_min,vmax=slip_max)
        #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=plt.cm.bone_r,vmin=slip_min,vmax=slip_max)
        #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=plt.cm.magma_r,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
        
        if afters==True:
            plt.scatter(along_strike_afters,depth_afters,marker='o',edgecolor='k',s=15,linewidth=0.5,facecolor='#CCCCCC')
        
        plt.ylabel('Depth (km)',fontsize=14)
        plt.xlabel('Along-strike distance (km)',fontsize=14)
        plt.xlim(xlims)
        plt.ylim(ylims)         
        plt.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='w',s=size,linewidth=2)
        for k in range(len(along_strike)):
            scale_slip=slip[k]/slip.max()
#            plt.quiver(along_strike[k],depth[k],rakess[k]/sqrt(rakess[k]**2+rakeds[k]**2),rakeds[k]/sqrt(rakess[k]**2+rakeds[k]**2),color='green',width=0.002,scale=50/scale_slip)
    #plt.annotate('B',xy=(-160,-4),fontsize=16,annotation_clip=False)
    #plt.annotate("B'",xy=(55,-4),fontsize=16,annotation_clip=False)
    cb.set_label('Slip (m)')
    #cb.set_label('CV')
    #cb.set_label('Std. Dev. (m)')
    plt.subplots_adjust(left=0.15, bottom=0.25, right=0.92, top=0.90, wspace=0, hspace=0)
    plt.grid(linestyle='dotted')
    plt.show()
    
    if histograms==True:
        plt.figure()
        plt.hist(depth_afters,50,orientation='horizontal',color='#FA8072')
        plt.ylim(ylims)

        plt.figure()
        plt.hist(along_strike_afters,50,color='#4682B4')
        plt.xlim(xlims)


def tile_resolution(rupt,resfile,nstrike,ndip,res_min,res_max,epicenter=0,epicenter_line=0):
    '''
    Plot subfault source time functions
    
    Parameters:
    epicenter is the coordinates, epcienter line is the down dip lien number where 
    the epcienter is
    '''
    
    from numpy import genfromtxt,unique,where,zeros,arange,pi,tile
    import matplotlib.pyplot as plt
    from obspy.core.util.geodetics import gps2DistAzimuth
    from matplotlib import cm
    
    f=genfromtxt(rupt)
    num=f[:,0]
    unum=unique(num)
    res=genfromtxt(resfile,usecols=2)*30
    res2=genfromtxt(u'/Users/dmelgar/Slip_inv/Napa_seis/analysis/resolution/seis_vonly_1winnpy.R',usecols=2)*5
    res=res+res2
    #Do same thing for centroid position
    loncent=-122.313
    latcent=38.26
    zcent=-4
    out=gps2DistAzimuth(epicenter[1],epicenter[0],latcent,loncent)
    xcent=-(out[0]/1000)
    #Done with centroid
    lon=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,1] #Only compute line at the epicenter depth
    lat=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,2]    
    depth=-f[:,3]
    depth=depth[0:len(unum)]
    along_strike=zeros(nstrike)
    for k in range(len(lat)):
        out=gps2DistAzimuth(epicenter[1],epicenter[0],lat[k],lon[k])
        if lat[k]<epicenter[1]: #It's to the south
            along_strike[k]=out[0]/1000
        else:
            along_strike[k]=-out[0]/1000
    #Now tile
    along_strike=tile(along_strike,ndip)
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1          
    #Plot
    plt.figure(num=None, figsize=(8, 4), dpi=80)
    plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=res,s=250,cmap=cm.gist_stern,vmin=res_min,vmax=res_max)
    cb=plt.colorbar()
    plt.ylabel('Depth (km)')
    plt.xlabel('Along-strike distance (km)')
    plt.xlim(along_strike.min()-1,along_strike.max()+1)
    plt.ylim(depth.min()-1,depth.max()+1)
    plt.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='#00FF00',s=350,linewidth=2)
    plt.scatter(xcent,zcent,marker='D',edgecolor='black',facecolor='',s=120,linewidth=2)
    cb.set_label('Resolution')
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.95, wspace=0, hspace=0)
    plt.show()


   
def tile_slip_movieframes(home,project_name,sliprate_path,nstrike,ndip,slip_min,slip_max,dt,geographic=False,epicenter=0,epicenter_line=0):
    '''
    Quick and dirty plot of a .rupt file
    epicenter is the coordinates, epcienter line is the down dip lien number where 
    the epcienter is
    '''
    
    from numpy import genfromtxt,zeros,tile,linspace,pi,cos,sin,ones,meshgrid,arange,r_,c_,sort
    import matplotlib.pyplot as plt
    from obspy.geodetics import gps2dist_azimuth
    from glob import glob
    from scipy.interpolate import griddata

    
    #Get sliprate files
    files=sort(glob(sliprate_path+'*.sliprate'))
    for kframe in range(len(files)):
        print(kframe)
        f=genfromtxt(files[kframe])
        slip=f[:,9]
        #Add aftershocks
        #afters=genfromtxt('/Users/dmelgar/Napa2014/hardebeck_afters.txt')
#        lonaf=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=4)
#        lataf=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=3)
#        zaf=-genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=5)
        #lonaf=afters[:,2]
        #lataf=afters[:,1]
        #zaf=-afters[:,3]
#        xaf=zeros(zaf.shape)
#        for k in range(len(lataf)):
#            out=gps2dist_azimuth(epicenter[1],epicenter[0],lataf[k],lonaf[k])
#            xaf[k]=(out[0]/1000)
#            if lataf[k]<epicenter[1]: #If it's tot he left of epcietner it's engative
#                xaf[k]=-xaf[k]
        #Done with afters
        #Do same thing for centroid position
        #loncent=-122.313
        #latcent=38.26
        #zcent=-4
        #out=gps2DistAzimuth(epicenter[1],epicenter[0],latcent,loncent)
        #xcent=-(out[0]/1000)
        #Done with centroid
        lon=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,1] #Only compute line at the epicenter depth
        lat=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,2]    
        depth=-f[:,3]
        along_strike=zeros(nstrike)
        for k in range(len(lat)):
            out=gps2dist_azimuth(epicenter[1],epicenter[0],lat[k],lon[k])
            if lat[k]>epicenter[1]: #It's to the south
                along_strike[k]=-out[0]/1000
            else:
                along_strike[k]=out[0]/1000
        #Now tile
        along_strike=tile(along_strike,ndip)
        #Get indices for plot
        istrike=zeros(nstrike*ndip)
        idip=zeros(nstrike*ndip)
        k=0
        for i in range(ndip):
            for j in range(nstrike):
                istrike[k]=nstrike-j
                idip[k]=ndip-i
                k+=1          
        #Make rupture velocity contours
        theta=linspace(0,2*pi,100)
        t=dt*kframe #Current time
        print("t="+str(t))
        r15=(2.2*t)*ones(100)
        x15=r15*cos(theta)
        y15=r15*sin(theta)-epicenter[2]
        r20=(2.0*t)*ones(100)
        x20=r20*cos(theta)
        y20=r20*sin(theta)-epicenter[2]
        r25=(2.5*t)*ones(100)
        x25=r25*cos(theta)
        y25=r25*sin(theta)-epicenter[2]
        r30=(3.0*t)*ones(100)
        x30=r30*cos(theta)
        y30=r30*sin(theta)-epicenter[2]
        #Plot
        plt.figure(num=None, figsize=(14, 3.5), dpi=80)
        #This plots slip as individual markers
        #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        #cb=plt.colorbar()
        #End single marker
        #This interpolates and plots the slip as a surface
        #First get limits of plot
        xlim=[along_strike.min()-0.5,along_strike.max()+0.5]
        ylim=[depth.min()-0.5,depth.max()+0.5]
        #Fix edges
        x1=along_strike[0:nstrike]
        y1=(depth[0]+0.5)*ones(x1.shape)
        z1=slip[0:nstrike]
        x2=along_strike[-nstrike:]
        y2=(depth[-1]-0.5)*ones(x2.shape)
        z2=slip[-nstrike:]
        x3=along_strike[arange(0,nstrike*ndip,nstrike)]+0.5
        y3=depth[arange(0,nstrike*ndip,nstrike)]
        z3=slip[arange(0,nstrike*ndip,nstrike)]
        x4=along_strike[arange(nstrike-1,nstrike*ndip,nstrike)]-0.5
        y4=depth[arange(nstrike-1,nstrike*ndip,nstrike)]
        z4=slip[arange(nstrike-1,nstrike*ndip,nstrike)]
        x5=along_strike[0]+0.5
        y5=depth[0]+0.5
        z5=slip[0]
        x6=along_strike[nstrike-1]-0.5
        y6=depth[nstrike-1]+0.5
        z6=slip[nstrike-1]
        x7=along_strike[-nstrike]+0.5
        y7=depth[-nstrike]-0.5
        z7=slip[nstrike-1]
        x8=along_strike[-1]-0.5
        y8=depth[-1]-0.5
        z8=slip[-1]
        along_strike=r_[along_strike,x1,x2,x3,x4,x5,x6,x7,x8]
        depth=r_[depth,y1,y2,y3,y4,y5,y6,y7,y8]
        slip=r_[slip,z1,z2,z3,z4,z5,z6,z7,z8]
        x=linspace(along_strike.min()-0.5,along_strike.max()+0.5,300)  #Adjsut for the width of the subfault
        y=linspace(depth.min()-0.49,depth.max()+0.49,100)   #Adjsut for height of subfault
        X, Y = meshgrid(x, y)
        sliprate=griddata(c_[along_strike,depth],slip,c_[X.ravel(),Y.ravel()],method='cubic')
        plt.scatter(X.ravel(),Y.ravel(),c=sliprate,marker='s',linewidth=0.0,edgecolor='',cmap=whitejet,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
#        plt.contourf(X,Y,sliprate,100,vmin=slip_min,vmax=slip_max,cmap=whitejet)
        plt.grid()
        #End interpolated 
        plt.ylabel('Depth (km)')
        plt.xlabel('Along-strike distance (km)')

        plt.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='#00FF00',s=350,linewidth=2)
        #plt.scatter(xcent,zcent,marker='D',edgecolor='black',facecolor='',s=120,linewidth=2)

##        plt.scatter(xaf,zaf,edgecolor='k',s=5)
        plt.plot(x15,y15,'--',c='grey')
#        plt.plot(x20,y20,'--',c='grey')
#        plt.plot(x25,y25,'--',c='grey')

        # grids
        major_ticks = arange(-100, 101, 10)
        minor_ticks = arange(-100, 101, 2)
        
        ax=plt.gca()
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        
        # And a corresponding grid
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)
        plt.xlim([-15,38])
        plt.ylim([-21,0])
    


        cb.set_label('Slip rate (m/s)')
        plot_name=home+project_name+'/plots/sliprate.'+str(kframe).rjust(4,'0')+'.pdf'
        plt.title('t = %.2fs' % (kframe*dt))
        plt.subplots_adjust(left=0.1, bottom=0.18, right=0.95, top=0.9, wspace=0, hspace=0)
        plt.annotate(s='NE',xy=(-14,-3),fontsize=18)
        plt.annotate(s='SW',xy=(35,-3),fontsize=18)
#        plt.axis('equal')
        plt.savefig(plot_name)     
#        ts=kframe*dt
#        if ts==1.0 or ts==2.0 or ts==3.0 or ts==4.0 or ts==5.0 or ts==6.0 or ts==7.0 or ts==8.0:
#            plot_name=home+project_name+'/plots/sliprate.'+str(kframe).rjust(4,'0')+'.pdf'
#            plt.savefig(plot_name) 
#            print('Saved PDF frame')
        plt.close("all")
    
        
            
def panel_tile_slip(home,project_name,sliprate_path,nstrike,ndip,slip_min,slip_max,nframes,geographic=False,epicenter=0,epicenter_line=0):
    '''
    Quick and dirty plot of a .rupt file
    epicenter is the coordinates, epcienter line is the down dip lien number where 
    the epcienter is
    
    nframes is [0,1,2,3]
    '''
    
    from numpy import genfromtxt,zeros,tile,linspace,pi,cos,sin,ones,meshgrid,arange,r_,reshape,where
    import matplotlib.pyplot as plt
    from obspy.core.util.geodetics import gps2DistAzimuth
    from glob import glob
    from matplotlib.mlab import griddata
    from matplotlib import colorbar
    
    #Get sliprate files
    files=[]
    for k in range(len(nframes)):
        frame=str(nframes[k]).rjust(4,'0')
        files.append(sliprate_path+frame+'.slip')
    #Now intialize figure
    fig, axarr = plt.subplots(5, 2,figsize=(8,4))
    axarr=axarr.transpose()
    axarr=reshape(axarr,(10,))
    #loop through frames
    for kframe in range(len(files)):
        ax=axarr[kframe] #current axis
        print(kframe)
        f=genfromtxt(files[kframe])
        slip=f[:,9]
        #Add aftershocks
        #afters=genfromtxt('/Users/dmelgar/Napa2014/hardebeck_afters.txt')
        lonaf=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=4)
        lataf=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=3)
        zaf=-genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=5)
        #lonaf=afters[:,2]
        #lataf=afters[:,1]
        #zaf=-afters[:,3]
        xaf=zeros(zaf.shape)
        for k in range(len(lataf)):
            out=gps2DistAzimuth(epicenter[1],epicenter[0],lataf[k],lonaf[k])
            xaf[k]=(out[0]/1000)
            if lataf[k]<epicenter[1]: #If it's tot he left of epcietner it's engative
                xaf[k]=-xaf[k]
        #Done with afters
        #Do same thing for centroid position
        #loncent=-122.313
        #latcent=38.26
        #zcent=-4
        #out=gps2DistAzimuth(epicenter[1],epicenter[0],latcent,loncent)
        #xcent=-(out[0]/1000)
        #Done with centroid
        lon=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,1] #Only compute line at the epicenter depth
        lat=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,2]    
        depth=-f[:,3]
        along_strike=zeros(nstrike)
        for k in range(len(lat)):
            out=gps2DistAzimuth(epicenter[1],epicenter[0],lat[k],lon[k])
            if lat[k]<epicenter[1]: #It's to the south
                along_strike[k]=-out[0]/1000
            else:
                along_strike[k]=out[0]/1000
        #Now tile
        along_strike=tile(along_strike,ndip)
        #Get indices for plot
        istrike=zeros(nstrike*ndip)
        idip=zeros(nstrike*ndip)
        k=0
        for i in range(ndip):
            for j in range(nstrike):
                istrike[k]=nstrike-j
                idip[k]=ndip-i
                k+=1          
        #Make rupture velocity contours
        theta=linspace(0,2*pi,100)
        t=1.0*(kframe+1) #Current time
        print("t="+str(t))
        r15=(1.5*t)*ones(100)
        x15=r15*cos(theta)
        y15=r15*sin(theta)-epicenter[2]
        r20=(2.0*t)*ones(100)
        x20=r20*cos(theta)
        y20=r20*sin(theta)-epicenter[2]
        r25=(2.5*t)*ones(100)
        x25=r25*cos(theta)
        y25=r25*sin(theta)-epicenter[2]
        r30=(3.0*t)*ones(100)
        x30=r30*cos(theta)
        y30=r30*sin(theta)-epicenter[2]
        #Plot
        
        #This plots slip as individual markers
        #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        #cb=plt.colorbar()
        #End single marker
        #This interpolates and plots the slip as a surface
        #First get limits of plot
        xlim=[along_strike.min()-0.5,along_strike.max()+0.5]
        ylim=[depth.min()-0.5,depth.max()+0.5]
        #Fix edges
        x1=along_strike[0:nstrike]
        y1=(depth[0]+0.5)*ones(x1.shape)
        z1=slip[0:nstrike]
        x2=along_strike[-nstrike:]
        y2=(depth[-1]-0.5)*ones(x2.shape)
        z2=slip[-nstrike:]
        x3=along_strike[arange(0,nstrike*ndip,nstrike)]+0.5
        y3=depth[arange(0,nstrike*ndip,nstrike)]
        z3=slip[arange(0,nstrike*ndip,nstrike)]
        x4=along_strike[arange(nstrike-1,nstrike*ndip,nstrike)]-0.5
        y4=depth[arange(nstrike-1,nstrike*ndip,nstrike)]
        z4=slip[arange(nstrike-1,nstrike*ndip,nstrike)]
        x5=along_strike[0]+0.5
        y5=depth[0]+0.5
        z5=slip[0]
        x6=along_strike[nstrike-1]-0.5
        y6=depth[nstrike-1]+0.5
        z6=slip[nstrike-1]
        x7=along_strike[-nstrike]+0.5
        y7=depth[-nstrike]-0.5
        z7=slip[nstrike-1]
        x8=along_strike[-1]-0.5
        y8=depth[-1]-0.5
        z8=slip[-1]
        along_strike=r_[along_strike,x1,x2,x3,x4,x5,x6,x7,x8]
        depth=r_[depth,y1,y2,y3,y4,y5,y6,y7,y8]
        slip=r_[slip,z1,z2,z3,z4,z5,z6,z7,z8]
        x=linspace(along_strike.min()-0.5,along_strike.max()+0.5,100)  #Adjsut for the width of the subfault
        y=linspace(depth.min()-0.49,depth.max()+0.49,100)   #Adjsut for height of subfault
        X, Y = meshgrid(x, y)
        sliprate=griddata(along_strike,depth,slip,X,Y,interp='linear')
        ax.scatter(along_strike,depth,marker='s',linewidth=0.0,edgecolor='',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        if kframe==5:
            i,j=where(sliprate.data>slip_max)
            sliprate.data[i,j]=slip_max
            cf=ax.contourf(X,Y,sliprate,100,vmin=slip_min,vmax=slip_max,cmap=whitejet)
        else:
            ax.contourf(X,Y,sliprate,100,vmin=slip_min,vmax=slip_max,cmap=whitejet)
        #ax.grid()
        #End interpolated 
        ax.set_xlim([-40.5,40.5])
        ax.set_ylim(ylim)
        ax.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='#00FF00',s=200,linewidth=1)
        #plt.scatter(xcent,zcent,marker='D',edgecolor='black',facecolor='',s=120,linewidth=2)

        ax.scatter(xaf,zaf,edgecolor='k',s=0.2)
        ax.plot(x15,y15,c='grey')
        ax.plot(x20,y20,c='grey')
        ax.plot(x25,y25,c='grey')
        
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
        anot=str(kframe)+'-'+str(kframe+1)+'s'
        ax.annotate(anot,xy=(25,-6))
        
        if kframe==4 or kframe==9:
            ax.xaxis.set_ticklabels(['','-40','','-20','','0','','20','','40'])
            ax.set_xlabel('Along strike distance (km)')
        if kframe==0 or kframe==1 or kframe==2 or kframe==3 or kframe==4:
            ax.yaxis.set_ticklabels(['','16','','','','','','','0'])
        if kframe==2:
            ax.set_ylabel('Down-dip distance (km)',rotation=90)
    
    plt.subplots_adjust(bottom=0.12,hspace=0.1,wspace=0.1,top=0.96)
    #Make stupid colorbar
    cax,kw = colorbar.make_axes([ax for ax in axarr.flat])
    cax.yaxis.set_ticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    plt.colorbar(cf, cax=cax,label='Slip(m)',ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6],**kw)
    plt.show() 
                    
                            
                       
def tile_moment(rupt,epicenter,nstrike,ndip,covfile,beta=0,vfast=0,vslow=0,shade=False,single_force=0):
    '''
    Tile plot of subfault source-time functions
    '''
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np
    from numpy import genfromtxt,unique,zeros,where,meshgrid,linspace,load,arange,expand_dims,squeeze,tile,r_
    from mudpy.forward import get_source_time_function,add2stf
    from mudpy.inverse import d2epi,ds2rot
    

    f=genfromtxt(rupt)
    num=f[:,0]
    nfault=nstrike*ndip
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    all=zeros(len(all_ss)*2)
    iss=2*arange(0,len(all)/2,1).astype(int)
    ids=2*arange(0,len(all)/2,1).astype(int)+1
    all[iss]=all_ss
    all[ids]=all_ds
    rot=ds2rot(expand_dims(all,1),beta)
    #Compute CI
    #Load covariances
    if covfile!=None:
        C=load(covfile)
        CIplus=squeeze(rot)+1*(C**0.5)
        CIminus=squeeze(rot)-1*(C**0.5)
        CIminus[CIminus<0]=0
        slipCIplus=(CIplus[iss]**2+CIplus[ids]**2)**0.5
        slipCIminus=(CIminus[iss]**2+CIminus[ids]**2)**0.5
    #Now parse for multiple rupture speeds
    unum=unique(num)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Get coordinates and compute distances
    source=f[0:len(unum),1:4]
    d=d2epi(epicenter,source)
    #Define velocity limits
    
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j-1
             idip[k]=i
             k+=1  
    #Define canvas
    fig, axarr = plt.subplots(ndip, nstrike)
    #Loop over subfaults
    Mmax=0
    Mout=[]
    for kfault in range(nfault):
        if kfault%10==0:
            print('... working on subfault '+str(kfault)+' of '+str(nfault))
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        if covfile !=None:
            slip_plus=slipCIplus[i]
            slip_minus=slipCIminus[i]
        #Get first source time function
        t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0],single_force=single_force)
        if covfile !=None:
            t1plus,M1plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_plus[0],single_force=single_force)
            t1minus,M1minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_minus[0],single_force=single_force)
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1],single_force=single_force)
            if covfile !=None:
                t2plus,M2plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_plus[kwin+1],single_force=single_force)
                t2minus,M2minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_minus[kwin+1],single_force=single_force)
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
            if covfile !=None:
                t1plus,M1plus=add2stf(t1plus,M1plus,t2plus,M2plus)
                t1minus,M1minus=add2stf(t1minus,M1minus,t2minus,M2minus)
        #Save M1 for output
        if kfault==0:
            Mout=expand_dims(M1,1).T
            tout=expand_dims(t1,1).T
        else:
            Mout=r_[Mout,expand_dims(M1,1).T]
            tout=r_[tout,expand_dims(t1,1).T]
        #Track maximum moment
        Mmax=max(Mmax,M1.max())
        #Done now plot them
        #get current axis
        ax=axarr[int(istrike[kfault])]
        #ax=axarr[int(idip[kfault]), int(istrike[kfault])] ###############################################
        if shade:
            #Make contourf
            Mc=linspace(0,0.98*max(M1),100)
            T,M=meshgrid(t1,Mc)
            i=where(T==0)[0]
            T[i]=0.01
            V=d[kfault]/T
            im=ax.contourf(T,M,V,100,vmin=vslow,vmax=vfast,cmap=cm.spectral)
            #Cover upper part
            ax.fill_between(t1,y1=M1,y2=1.01*M1.max(),color='white')
            #Plot confidence intervals
            if covfile !=None:
                ax.fill_between(t1,M1minus,M1plus,facecolor='grey',alpha=0.4)
                ax.plot(t1,M1plus,color='black')
                ax.plot(t1,M1minus,color='white',lw=2)
        #Plot curve
##########################################
        #plt.figure(kfault)
        #plt.title(f'{kfault}')
        #plt.plot(t1,M1)
#######################################
        ax.plot(t1, M1,color='k')
        
        ax.grid()
        ax.set_xlim([t1[0],t1[-1]])
        ax.xaxis.set_ticks(linspace(t1[0],t1[-1],5))
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    #Go back and rescale all subplots by maximum moment
    for k in range(ndip):
        for k2 in range(nstrike):
            ax=axarr[k]
            #ax=axarr[k,k2]###########################################
            ax.set_ylim([0,Mmax])
    #Fix subplot arrangement
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.9, top=0.98, wspace=0, hspace=0)
    #Add colorbar
    if shade:
        T2=T.copy()
        M2=M.copy()
        M2[:,:]=0
        vfill=linspace(1.0,3.0,10000)
        V2=tile(vfill[::-1],(100,1))
        ax=axarr[0,0]
        im2=ax.contourf(T2,M2,V2,100,vmin=vslow,vmax=vfast,cmap=cm.spectral)
        cbar_ax = fig.add_axes([0.91, 0.15, 0.01, 0.7])
        cb=fig.colorbar(im2, cax=cbar_ax)
        cb.set_label('Reference rupture velocity (km/s)')
    print('Maximum moment was '+str(Mmax)+'N-m')
    return tout,Mout


def source_time_function(rupt,epicenter,plot=True,xlim=None,ylim=None,normalize=True,single_force=0):
    '''
    Plot source time function of complete rupture
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,unique,log10,where,floor
    from mudpy.forward import get_source_time_function,add2stf
    
    f=genfromtxt(rupt)
    num=f[:,0]
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    nfault=len(unum)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print('... working on subfault '+str(kfault)+' of '+str(nfault))

        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        if kfault==0:#Get first source time function
            t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0],single_force=single_force)
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1],single_force=single_force)
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
    #Get power
    exp=floor(log10(M1.max()))
    if normalize == True:
        M1=M1/(10**exp)
    else:
        pass
    
    if plot==True:
        plt.figure()
        plt.fill(t1,M1,'b',alpha=0.5)
        plt.plot(t1,M1,color='k')
        plt.grid()
        plt.xlabel('Time(s)')
        if single_force == 1:
            plt.ylabel(f'Force Rate {exp} (Dyne or Newton/s)')
        else:
            plt.ylabel('Moment Rate ('+r'$\times 10^{'+str(int(exp))+r'}$Nm/s)')
        plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.7, wspace=0, hspace=0)
        if xlim!=None:
            plt.xlim(xlim)
        if ylim!=None:
            plt.ylim(ylim)
    return t1,M1
    

def source_time_function_FQ(rupt,epicenter,plot=False):
    '''
    Get source time function for Fakequake (fix bug in rise_time=0 in subfaults)
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,unique,log10,where,floor
    from mudpy.forward import get_source_time_function,add2stf
    
    f=genfromtxt(rupt)
    num=f[:,0]
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    nfault=len(unum)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Loop over subfaults
    t1=0;M1=0 #reset t1 and M1
    for kfault in range(nfault):
        if kfault%10==0:
            print('... working on subfault '+str(kfault)+' of '+str(nfault))
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        if type(M1)==int and rise_time[kfault]!=0: #Get first source time function
            t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        elif type(M1)!=int and rise_time[kfault]!=0:
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
            t1,M1=add2stf(t1,M1,t2,M2)
    #Get power
    exp=floor(log10(M1.max()))
    M2=M1.copy()
    M1=M1/(10**exp)
    if plot==True:
        plt.figure()
        plt.fill(t1,M1,'b',alpha=0.5)
        plt.plot(t1,M1,color='k')
        plt.grid()
        plt.xlabel('Time(s)')
        plt.ylabel('Moment Rate ('+r'$\times 10^{'+str(int(exp))+r'}$Nm/s)')
        plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.7, wspace=0, hspace=0)
    return t1,M1,M2



def geographic_STFs(rupt,epicenter,nstrike,ndip,tscale=100,Mscale=1,figsize=(8,10),tout=[],Mout=[]):
    '''
    Plot STFs in their geographic locations
    '''
    
    from mudpy import analysis
    from numpy import genfromtxt
    from matplotlib import pyplot as plt
    #Read source file; and determine subfault coordinate
    fault=genfromtxt(rupt)
    N=ndip*nstrike
    fault=fault[0:N,:]
    lon=fault[:,1]
    lat=fault[:,2]
    if tout==[] or Mout==[]:
        tout,Mout=analysis.subfault_STFs(rupt,epicenter,223,1)
    #Determine peak moment rate
    Mmax=Mout.max()
    Mout=Mout/Mmax
    #determine plot limits
    lon_extra=(lon.max()-lon.min())
    lat_extra=(lat.max()-lat.min())
    #Start plot
    plt.figure(figsize=figsize)
    plt.scatter(lon,lat)
    plt.scatter(epicenter[0],epicenter[1],marker='*',c='r',s=200)
    #
    for k in range(len(lon)):
        t1=tout[k,0]
        tplot=((tout[k,:]-t1)/tscale)+lon[k]
        Mplot=(Mout[k,:]/Mscale)+lat[k]
        plt.plot(tplot,Mplot,'k')
    plt.title(rupt.split('/')[-1])
    plt.show()
    
    


def Rm(G,lambda_spatial,lambda_temporal,Ls,Lt,bounds,nstrike,ndip,maxR=0.2):
    '''
    Plot model resolution matrix
    '''
    
    from numpy import diag,zeros,arange
    import matplotlib.pyplot as plt
    from numpy.linalg import inv
    
    #Compute model resolution matrix
    Gs=G.transpose().dot(G)+(lambda_spatial**2)*Ls.transpose().dot(Ls)+(lambda_temporal**2)*Lt.transpose().dot(Lt)
    R=inv(Gs).dot(G.transpose()).dot(G)
    #Keep diagonals only
    r=diag(R)
    ids=arange(1,len(r),2)
    r=r[ids]
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1           
    #Plot
    plt.figure()
    plt.scatter(istrike,idip,marker='o',c=r,s=250,cmap=plt.cm.PuBuGn,vmax=maxR)
    plt.ylabel('Along-dip index')
    plt.xlabel('Along-strike index')
    cb=plt.colorbar()
    cb.set_label('Diagonal Value of R')
    plt.axis('equal')
    plt.xlim(istrike.min()-1,istrike.max()+1)
    plt.ylim(idip.min()-1,idip.max()+1)
    plt.grid()
    plt.title('Model Resolution')
    plt.show()
    return R

def tslice(rupt,out,dt,cumul):
    '''
    Quick and dirty plot of a .rupt file
    '''
    
    from numpy import genfromtxt,unique,where,zeros,arange,intersect1d,trapz
    import matplotlib.pyplot as plt
    
    delta_t=0.05
    f=genfromtxt(rupt)
    trupt=f[:,12]
    trise=f[:,7]
    all_ss=f[:,8]
    all_ds=f[:,9]
    num=f[:,0]
    #Get other parameters
    #lon=f[0:len(unum),1]
    #lat=f[0:len(unum),2]
    #strike=f[0:len(unum),4]
    #Decide on time vector
    tslice=arange(0,trupt.max()+dt,dt)
    #Determine triangle height at all subfaults
    hss=2*all_ss/trise
    hds=2*all_ds/trise
    #Cumulative
    ss_cumul=zeros(len(f))
    ds_cumul=zeros(len(f))
    #Determine time series fo each triangle
    t=arange(0,trupt.max()+trise[0],delta_t)
    for kslice in range(len(tslice-2)):
        print(str(kslice)+'/'+str(len(tslice)-1))
        #And initalize slice vectors
        ss_slice=zeros(len(f))
        ds_slice=zeros(len(f))
        for kfault in range(len(f)):
            yss=zeros(t.shape)
            yds=zeros(t.shape)
            #Up going
            i1=where(t>=trupt[kfault])[0]
            i2=where(t<=(trupt[kfault]+trise[0]/2))[0] #Ascending triangle
            i=intersect1d(i1,i2)
            yss[i]=(2*hss[kfault]/trise[0])*t[i]-(2*hss[kfault]*trupt[kfault]/trise[0])
            yds[i]=(2*hds[kfault]/trise[0])*t[i]-(2*hds[kfault]*trupt[kfault]/trise[0])
            #Down going
            i1=where(t>(trupt[kfault]+trise[0]/2))[0]
            i2=where(t<=(trupt[kfault]+trise[0]))[0] #Ascending triangle
            i=intersect1d(i1,i2)
            yss[i]=(-2*hss[kfault]/trise[0])*t[i]+(2*hss[kfault]/trise[0])*(trupt[kfault]+trise[0])
            yds[i]=(-2*hds[kfault]/trise[0])*t[i]+(2*hds[kfault]/trise[0])*(trupt[kfault]+trise[0])
            #Now integrate slip at pertinent time interval
            i1=where(t>=tslice[kslice])[0]
            i2=where(t<=tslice[kslice+1])[0]
            i=intersect1d(i1,i2)
            ss_slice[kfault]=trapz(yss[i],t[i])
            ds_slice[kfault]=trapz(yds[i],t[i])
        #Combine into single model for that time slice
        ss_cumul=ss_cumul+ss_slice
        ds_cumul=ds_cumul+ds_slice
        unum=unique(num)
        lon=f[0:len(unum),1]
        lat=f[0:len(unum),2]
        strike=f[0:len(unum),4]
        ss=zeros(len(unum))
        ds=zeros(len(unum))
        for k in range(len(unum)):
            if cumul==0:
                i=where(unum[k]==num)
                ss[k]=ss_slice[i].sum()
                ds[k]=ds_slice[i].sum()    
            else:
                i=where(unum[k]==num)
                ss[k]=ss_cumul[i].sum()
                ds[k]=ds_cumul[i].sum()        
        slip=(ss**2+ds**2)**0.5
        #Plot
        #Get projection of rake vector
        x,y=slip2geo(ss,ds,strike)
        #Plot
        plt.figure()
        plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=plt.cm.gnuplot2_r,vmin=0,vmax=35)
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        cb=plt.colorbar()
        plt.quiver(lon,lat,x,y,color='green',width=0.0013)
        plt.grid()
        if cumul==0:
            cb.set_label('Slip (m)')
            plt.title('t = '+str(tslice[kslice])+'s to '+str(tslice[kslice+1])+'s') 
            plt.savefig(out+str(kslice).rjust(4,'0')+'.kin_slice.png')
        else:
            cb.set_label('Cumulative Slip (m)')
            plt.title('t = '+str(tslice[kslice+1])+'s')
            plt.savefig(out+str(kslice).rjust(4,'0')+'.kin_cumulative.png')
        plt.close("all")
    

def plot_data(home,project_name,gflist,vord,decimate,lowpass,t_lim,sort,scale,k_or_g):
    '''
    Plot ( vs real data
    
    gflist: The GF control fiel that decides what to plot/not plot
    datapath
    '''
    from obspy import read
    from numpy import genfromtxt,where,argsort
    import matplotlib.pyplot as plt
    import matplotlib
    from mudpy.green import stdecimate 
    from mudpy.forward import lowpass as lfilter
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='S')
    lon=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/data/waveforms/'
    if vord.lower()=='d':
        kgf=0 #disp
        if k_or_g.lower()=='kal':
            datasuffix='kdisp'
        else:
            datasuffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        datasuffix='vel'
    elif vord.lower()=='a':
        kgf=1 #disp
        datasuffix='acc'
    #Decide on sorting
    i=where(gf[:,kgf]==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    nsta=len(i)
    fig, axarr = plt.subplots(nsta, 3)  
    for k in range(len(i)):
        n=read(datapath+sta[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+sta[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+sta[i[k]]+'.'+datasuffix+'.u')
        if lowpass!=None:
            fsample=1./e[0].stats.delta
            e[0].data=lfilter(e[0].data,lowpass,fsample,2)
            n[0].data=lfilter(n[0].data,lowpass,fsample,2)
            u[0].data=lfilter(u[0].data,lowpass,fsample,2)
        if decimate!=None:
            n[0]=stdecimate(n[0],decimate)
            e[0]=stdecimate(e[0],decimate)
            u[0]=stdecimate(u[0],decimate)
        if scale!=None:
            n[0].data=n[0].data/scale
            e[0].data=e[0].data/scale
            u[0].data=u[0].data/scale
        #Make plot
        if nsta>1:
            axn=axarr[k,0]
            axe=axarr[k,1]
            axu=axarr[k,2]
        else:
            axn=axarr[0]
            axe=axarr[1]
            axu=axarr[2]
        axn.plot(n[0].times(),n[0].data,'k')
        axn.grid(which='both')
        axe.plot(e[0].times(),e[0].data,'k')
        axe.grid(which='both')
        axu.plot(u[0].times(),u[0].data,'k')
        axu.grid(which='both')
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axe.set_xlim(t_lim)
        axn.set_xlim(t_lim)
        axu.set_xlim(t_lim)
        axn.yaxis.set_ticklabels([])
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axn.yaxis.grid(False)
        axe.yaxis.grid(False)
        axu.yaxis.grid(False)
        axn.yaxis.set_ticks([])
        axe.yaxis.set_ticks([])
        axu.yaxis.set_ticks([])
        
        #Annotations
        trange=t_lim[1]-t_lim[0]
        sign=1.
        if abs(min(n[0].data))>max(n[0].data):
            sign=-1. 
        nmax='%.3f' % (sign*max(abs(n[0].data)))
        sign=1.
        nlims=axn.get_ylim()
        nrange=nlims[1]-nlims[0]
        
        if abs(min(e[0].data))>max(e[0].data):
            sign=-1.         
        emax='%.3f' % (sign*max(abs(e[0].data)))
        sign=1.
        elims=axe.get_ylim()
        erange=elims[1]-elims[0]
        
        if abs(min(u[0].data))>max(u[0].data):
            sign=-1. 
        umax='%.3f' % (sign*max(abs(u[0].data)))
        sign=1.
        ulims=axu.get_ylim()
        urange=ulims[1]-ulims[0]
        
        #axn.annotate(nmax,xy=(t_lim[1]-0.3*trange,nlims[0]+0.02*nrange),fontsize=12)
        #axe.annotate(emax,xy=(t_lim[1]-0.3*trange,elims[0]+0.02*erange),fontsize=12)
        #axu.annotate(umax,xy=(t_lim[1]-0.3*trange,ulims[0]+0.02*urange),fontsize=12)
        #axn.annotate(nsmax,xy=(t_lim[1]-0.3*trange,nlims[0]+0.7*nrange),fontsize=12,color='red')
        #axe.annotate(esmax,xy=(t_lim[1]-0.3*trange,elims[0]+0.7*erange),fontsize=12,color='red')
        #axu.annotate(usmax,xy=(t_lim[1]-0.3*trange,ulims[0]+0.7*urange),fontsize=12,color='red')
        axn.annotate(nmax,xy=(t_lim[1]-0.25*trange,nlims[0]+0.02*nrange),fontsize=12)
        axe.annotate(emax,xy=(t_lim[1]-0.25*trange,elims[0]+0.02*erange),fontsize=12)
        axu.annotate(umax,xy=(t_lim[1]-0.25*trange,ulims[0]+0.02*urange),fontsize=12)
        #Station name
        axn.set_ylabel(sta[i[k]],rotation=90)
        if k==0:
            if vord.lower()=='d':
                axn.set_title('North (m)')
                axe.set_title('East (m)')
                axu.set_title('Up (m)')
            elif vord.lower()=='v':
                axn.set_title('North (m/s)')
                axe.set_title('East (m/s)')
                axu.set_title('Up (m/s)')
            else:
                axn.set_title(r'North (m/s$^2$)')
                axe.set_title('East (m/s$^2$)')
                axu.set_title('Up (m/s$^2$)')
        if k!=len(i)-1 or len(i)==1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
            xtick=axn.xaxis.get_majorticklocs()
            #ix=[1,3,5]
            #ix=[2,4,6]
            ix=[1,3,5,7]
            ix=0
            xtick=xtick[ix]
            #xticklabel=['','50','','150','','250',''] #Tohoku
            #xticklabel=['0','','20','','40','','60'] #Napa preferred
            #xticklabel=['','10','','30','','50','','70'] #Napa preferred
            xticklabel=['','','40','','80','','120','','160',''] #Iquique preferred
        if k==len(i)-1 and nsta>1: #Last plot
            axe.set_xlabel('Time (s)')
            axn.xaxis.set_ticklabels(xticklabel)
            axe.xaxis.set_ticklabels(xticklabel)
            axu.xaxis.set_ticklabels(xticklabel)
            #axn.xaxis.set_ticks(xtick)
            #axe.xaxis.set_ticks(xtick)
            #axu.xaxis.set_ticks(xtick)
    #plt.subplots_adjust(left=0.2, bottom=0.05, right=0.8, top=0.95, wspace=0, hspace=0)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.8, top=0.85, wspace=0, hspace=0)
      
    
    
    
    
def synthetics(home,project_name,run_name,run_number,gflist,vord,decimate,lowpass,t_lim,
               sort,scale,k_or_g,uncert=False,waveforms_as_accel=False,units='m',uncerth=0.01,uncertv=0.03,
               tick_frequency=10,spoof_vel_as_disp=False,return_vectors=False,same_channel_amplitude=True,
               show=False):
    '''
    Plot synthetics vs real data
    
    gflist: The GF control fiel that decides what to plot/not plot
    datapath
    '''
    from obspy import read
    from numpy import genfromtxt,where,argsort,zeros,array,r_
    import matplotlib.pyplot as plt
    import matplotlib
    from mudpy.green import stdecimate 
    from mudpy.forward import lowpass as lfilter
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='U')
    lon=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    if vord.lower()=='d':
        kgf=0 #disp
        if k_or_g.lower()=='kal':
            datasuffix='kdisp'
        else:
            datasuffix='disp'
        synthsuffix='disp'
        if spoof_vel_as_disp==True: #Your an SM as dispalcememnt but was save as vel.
            synthsuffix='vel'
    elif vord.lower()=='v':
        kgf=1 #vel
        if waveforms_as_accel==False:
            datasuffix='vel'
        else:
            datasuffix='accel'
        synthsuffix='vel'
    #Decide on sorting
    i=where(gf[:,kgf]==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    nsta=len(i)
    fig, axarr = plt.subplots(nsta, 3)  
    dvector=array([])
    dsvector=array([])
    for k in range(len(i)):
        n=read(datapath+sta[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+sta[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+sta[i[k]]+'.'+datasuffix+'.u')
        ns=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.n.sac')
        es=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.e.sac')
        us=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.u.sac')
        
        if units=='cm':
            n[0].data*=100
            e[0].data*=100
            u[0].data*=100
            ns[0].data*=100
            es[0].data*=100
            us[0].data*=100
           
        
        if lowpass!=None:
            print('Lowpassing')
            fsample=1./e[0].stats.delta
            e[0].data=lfilter(e[0].data,lowpass,fsample,2)
            n[0].data=lfilter(n[0].data,lowpass,fsample,2)
            u[0].data=lfilter(u[0].data,lowpass,fsample,2)
            es[0].data=lfilter(es[0].data,lowpass,fsample,2)
            ns[0].data=lfilter(ns[0].data,lowpass,fsample,2)
            us[0].data=lfilter(us[0].data,lowpass,fsample,2)
        if decimate!=None:
            n[0]=stdecimate(n[0],decimate)
            e[0]=stdecimate(e[0],decimate)
            u[0]=stdecimate(u[0],decimate)

        if scale!=None:
            n[0].data=n[0].data/scale
            ns[0].data=ns[0].data/scale
            e[0].data=e[0].data/scale
            es[0].data=es[0].data/scale
            u[0].data=u[0].data/scale
            us[0].data=us[0].data/scale
        #Make plot
        if nsta>1:
            axn=axarr[k,0]
            axe=axarr[k,1]
            axu=axarr[k,2]
        else:
            axn=axarr[0]
            axe=axarr[1]
            axu=axarr[2]
            
        #add to output vector
        dvector = r_[dvector,n[0].data,e[0].data,u[0].data]
        dsvector = r_[dsvector,ns[0].data,es[0].data,us[0].data]
        
        axn.plot(n[0].times(),n[0].data,'k',ns[0].times(),ns[0].data,'r')
        axn.grid(which='both')
        axe.plot(e[0].times(),e[0].data,'k',es[0].times(),es[0].data,'r')
        axe.grid(which='both')
        axu.plot(u[0].times(),u[0].data,'k',us[0].times(),us[0].data,'r')
        if uncert==True:
            axn.fill_between(n[0].times(),n[0].data-uncerth,n[0].data+uncerth,alpha=0.2,color='k')
            axe.fill_between(e[0].times(),e[0].data-uncerth,e[0].data+uncerth,alpha=0.2,color='k')
            axu.fill_between(u[0].times(),u[0].data-uncertv,u[0].data+uncertv,alpha=0.2,color='k')
        axu.grid(which='both')
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axe.set_xlim(t_lim)
        axn.set_xlim(t_lim)
        axu.set_xlim(t_lim)
        axn.yaxis.set_ticklabels([])
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axn.yaxis.grid(False)
        axe.yaxis.grid(False)
        axu.yaxis.grid(False)
        axn.yaxis.set_ticks([])
        axe.yaxis.set_ticks([])
        axu.yaxis.set_ticks([])
        
        if same_channel_amplitude == True:
            ymax = max(abs(n[0].data.max()), abs(ns[0].data.max()), 
                       abs(e[0].data.max()), abs(es[0].data.max()), 
                       abs(u[0].data.max()), abs(us[0].data.max()))
            ymin = -ymax
            axn.set_ylim([ymin,ymax])
            axe.set_ylim([ymin,ymax])
            axu.set_ylim([ymin,ymax])

        #Annotations
        trange=t_lim[1]-t_lim[0]
        sign=1.
        if abs(min(n[0].data))>max(n[0].data):
            sign=-1. 
        nmax='%.2f' % (sign*max(abs(n[0].data)))
        sign=1.
        if abs(min(ns[0].data))>max(ns[0].data):
            sign=-1. 
        nsmax='%.2f' % (sign*max(abs(ns[0].data)))
        sign=1.
        nlims=axn.get_ylim()
        nrange=nlims[1]-nlims[0]
        
        if abs(min(e[0].data))>max(e[0].data):
            sign=-1.         
        emax='%.2f' % (sign*max(abs(e[0].data)))
        sign=1.
        if abs(min(es[0].data))>max(es[0].data):
            sign=-1. 
        esmax='%.2f' % (sign*max(abs(es[0].data)))
        sign=1.
        elims=axe.get_ylim()
        erange=elims[1]-elims[0]
        
        if abs(min(u[0].data))>max(u[0].data):
            sign=-1. 
        umax='%.2f' % (sign*max(abs(u[0].data)))
        sign=1.
        if abs(min(us[0].data))>max(us[0].data):
            sign=-1 
        usmax='%.2f' % (sign*max(abs(us[0].data)))
        sign=1.
        ulims=axu.get_ylim()
        urange=ulims[1]-ulims[0]
        
        #axn.annotate(nmax,xy=(t_lim[1]-0.3*trange,nlims[0]+0.02*nrange),fontsize=12)
        #axe.annotate(emax,xy=(t_lim[1]-0.3*trange,elims[0]+0.02*erange),fontsize=12)
        #axu.annotate(umax,xy=(t_lim[1]-0.3*trange,ulims[0]+0.02*urange),fontsize=12)
        #axn.annotate(nsmax,xy=(t_lim[1]-0.3*trange,nlims[0]+0.7*nrange),fontsize=12,color='red')
        #axe.annotate(esmax,xy=(t_lim[1]-0.3*trange,elims[0]+0.7*erange),fontsize=12,color='red')
        #axu.annotate(usmax,xy=(t_lim[1]-0.3*trange,ulims[0]+0.7*urange),fontsize=12,color='red')
        axn.annotate(nmax,xy=(t_lim[1]-0.25*trange,nlims[0]+0.02*nrange),fontsize=11)
        axe.annotate(emax,xy=(t_lim[1]-0.25*trange,elims[0]+0.02*erange),fontsize=11)
        axu.annotate(umax,xy=(t_lim[1]-0.25*trange,ulims[0]+0.02*urange),fontsize=11)
        axn.annotate(nsmax,xy=(t_lim[1]-0.25*trange,nlims[0]+0.7*nrange),fontsize=11,color='red')
        axe.annotate(esmax,xy=(t_lim[1]-0.25*trange,elims[0]+0.7*erange),fontsize=11,color='red')
        axu.annotate(usmax,xy=(t_lim[1]-0.25*trange,ulims[0]+0.7*urange),fontsize=11,color='red')
        #Station name
        #axn.set_ylabel(sta[i[k]],rotation=90)
        axn.set_ylabel(sta[i[k]],rotation=0,labelpad=20)
        
        #tick frequency
        axn.xaxis.set_minor_locator(MultipleLocator(tick_frequency))
        axe.xaxis.set_minor_locator(MultipleLocator(tick_frequency))
        axu.xaxis.set_minor_locator(MultipleLocator(tick_frequency))
        axn.xaxis.set_major_locator(MultipleLocator(tick_frequency))
        axe.xaxis.set_major_locator(MultipleLocator(tick_frequency))
        axu.xaxis.set_major_locator(MultipleLocator(tick_frequency))
        
        if k==0:
            if vord.lower()=='d':
                if units=='cm':
                    axn.set_title('North (cm)')
                    axe.set_title('East (cm)')
                    axu.set_title('Up (cm)')
                else:
                    axn.set_title('North (m)')
                    axe.set_title('East (m)')
                    axu.set_title('Up (m)')
            else:
                if units=='cm':
                    axn.set_title('North (cm/s)')
                    axe.set_title('East (cm/s)')
                    axu.set_title('Up (cm/s)')  
                else:
                    axn.set_title('North (m/s)')
                    axe.set_title('East (m/s)')
                    axu.set_title('Up (m/s)')
        if k!=len(i)-1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
#            xtick=axn.xaxis.get_majorticklocs()
#            #ix=[1,3,5]
#            #ix=[2,4,6]
#            ix=[1,3,5,7]
#            ix=0
#            xtick=xtick[ix]
            #xticklabel=['','50','','150','','250',''] #Tohoku
            #xticklabel=['0','','20','','40','','60'] #Napa preferred
            #xticklabel=['','10','','30','','50','','70'] #Napa preferred
            #xticklabel=['','100','','200','','300'] #Maule preferred
            #xticklabel=['','','40','','80','','120','','160',''] #Iquique preferred
            #xticklabel=['','10','','30','','50',''] #Nepal preferred
            #xticklabel=['','5','','15','','25','','35',''] #Lefkada prevferred
            #xticklabel=['','5','','15','','25',''] #Amatrice preferred
            #xticklabel=['','10','','30','','50','','70',''] #Melinka preferred
            #xticklabel=['','10','','30','','50','','70',''] #Tehuantepec preferred
#            xticklabel=['','20','','40','','60',''] #Tehuantepec preferred
#            xticklabel=axn.xaxis.get_ticklabels()
        if k==len(i)-1 and nsta>1: #Last plot
            axe.set_xlabel('Seconds since OT')
            xticklabel=axn.xaxis.get_ticklabels()
            print(xticklabel)
#            xticklabel=['','20','','40','','60','']
#            axn.set_xticks([20,40,60,80,100])
#            axe.set_xticks([20,40,60,80,100])
#            axu.set_xticks([20,40,60,80,100])
            # axn.set_xticks([15,30,45])
            # axe.set_xticks([15,30,45])
            # axu.set_xticks([15,30,45])
#            axn.xaxis.set_ticklabels(xticklabel)
#            axe.xaxis.set_ticklabels(xticklabel)
#            axu.xaxis.set_ticklabels(xticklabel)
            #axn.xaxis.set_ticks(xtick)
            #axe.xaxis.set_ticks(xtick)
            #axu.xaxis.set_ticks(xtick)
    #plt.subplots_adjust(left=0.2, bottom=0.05, right=0.8, top=0.95, wspace=0, hspace=0)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.8, top=0.85, wspace=0, hspace=0)
    if show:
        plt.show()
    
    if return_vectors == True:
        
        return dvector,dsvector

def static_synthetics(home,project_name,run_name,run_number,gflist,qscale,xl=None,yl=None):
    '''
    Plot synthetics vs real data
    
    gflist: The GF control fiel that decides what to plot/not plot
    datapath
    sscale: scales the synthetics if some weight has been applied
    qscale: scale fo the quiver plot arrows
    '''
    from numpy import genfromtxt,where,zeros
    import matplotlib.pyplot as plt
    import matplotlib
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='U')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[3],dtype='f')
    datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    n=zeros(len(i))
    e=zeros(len(i))
    u=zeros(len(i))
    ns=zeros(len(i))
    es=zeros(len(i))
    us=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(datapath+str(sta[i[k]])+'.neu')
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        n[k]=neu[0] ; e[k]=neu[1] ; u[k]=neu[2]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.static.neu')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        ns[k]=neus[0] ; es[k]=neus[1] ; us[k]=neus[2]
    #Make plot   
    lonadjust=(lon.max()-lon.min())/10
    latadjust=(lat.max()-lat.min())/10
    plt.figure()
    plt.subplot(121)
    plt.quiver(lon,lat,e,n,color='k',scale=qscale)
    plt.quiver(lon,lat,es,ns,color='r',scale=qscale)
    plt.grid()
    plt.title('Horizontals')
    #for k in range(len(i)):
    #    plt.annotate(sta[i[k]],xy=(lon[k],lat[k]))
    #plt.xlim([lon.min()-lonadjust,lon.max()+lonadjust])
    #plt.ylim([lat.min()-latadjust,lat.max()+latadjust])
    plt.xlim(xl)
    plt.ylim(yl)
    plt.subplot(122)
    plt.quiver(lon,lat,zeros(len(u)),u,color='k',scale=qscale)
    plt.quiver(lon,lat,zeros(len(us)),us,color='r',scale=qscale)
    plt.grid()
    plt.title('Verticals')
    #for k in range(len(i)):
    #    plt.annotate(sta[i[k]],xy=(lon[k],lat[k]))
    #plt.xlim([lon.min()-lonadjust,lon.max()+lonadjust])
    #plt.ylim([lat.min()-latadjust,lat.max()+latadjust])
    plt.xlim(xl)
    plt.ylim(yl)
    #plt.legend('Data','Synth')
    plt.suptitle('Statics for run '+project_name+': '+run_name+'.'+run_number)
    
def insar_residual(home,project_name,run_name,run_number,gflist,zlims):
    '''
    Plot insar_residual
    
    gflist: The GF control file that decides what to plot/not plot
    datapath
    '''
    from numpy import genfromtxt,where,zeros,sqrt,c_,savetxt
    import matplotlib.pyplot as plt
    import matplotlib
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='U')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[7],dtype='f')
    gf_datapath=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[12],dtype='U')
    #datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    gf_datapath=gf_datapath[i]
    los_data=zeros(len(i))
    los_synth=zeros(len(i))
    for k in range(len(i)):
        #neu=genfromtxt(datapath+sta[i[k]]+'.los')
        neu=genfromtxt(gf_datapath[k])
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        los_data[k]=neu[0]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.los')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        los_synth[k]=neus[0]
    #Make plot
    out=c_[lon,lat,los_data,los_synth,los_data-los_synth]
    savetxt(home+project_name+'/analysis/'+run_name+'.'+run_number+'.insar.res',out,fmt='%.6f\t%.6f\t%8.5f\t%8.5f\t%8.5f',header='lon,lat,los_data(m),los_synthetic(m),data-synthetic(m)')
    plt.figure()
    plt.scatter(lon,lat,c=los_data-los_synth,cmap=matplotlib.cm.seismic,vmin=zlims[0],vmax=zlims[1],s=5)
    plt.title('LOS data - LOS predicted (m)')
    plt.colorbar()
    plt.grid()
    
    
def insar_results(home,project_name,run_name,run_number,gflist,zlims,cmap,figsize=(8,5),title=None,interpolate=True,method='linear',npts=100,
                  show=False):
    '''
    Plot insar observed in one panel and insar modeled in the other
    
    gflist: The GF control file that decides what to plot/not plot
    datapath
    '''
    from numpy import genfromtxt,where,zeros,sqrt,c_,savetxt,linspace,meshgrid
    import matplotlib.pyplot as plt
    import matplotlib
    from scipy.interpolate import griddata
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='U')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[7],dtype='f')
    gf_datapath=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[12],dtype='U')
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    gf_datapath=gf_datapath[i]
    los_data=zeros(len(i))
    los_synth=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(gf_datapath[k])
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        los_data[k]=neu[0]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.los')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        los_synth[k]=neus[0]
    #Make plot
    out=c_[lon,lat,los_data,los_synth,los_data-los_synth]
    savetxt(home+project_name+'/analysis/'+run_name+'.'+run_number+'.insar.res',out,fmt='%.6f\t%.6f\t%8.5f\t%8.5f\t%8.5f',header='lon,lat,los_data(m),los_synthetic(m),data-synthetic(m)')
    
    plt.figure()
    
    #Interpolate to regular grid?
    if interpolate==True:
        x=linspace(lon.min(),lon.max(),npts)
        y=linspace(lat.min(),lat.max(),npts)
        X,Y=meshgrid(x,y)
        los_data_interp=griddata((lon,lat),los_data,(X,Y),method=method,fill_value=0)
        los_synth_interp=griddata((lon,lat),los_synth,(X,Y),method=method,fill_value=0)
        
        ax=plt.subplot(211)
        ax.tick_params(labelbottom='off') 
        plt.pcolormesh(X,Y,los_data_interp,cmap=cmap,vmin=zlims[0],vmax=zlims[1])
        plt.title('LOS observed (m)')
        plt.colorbar()
        plt.grid()
        plt.ylabel('Latitude')
        #replt.axis('equal')
        plt.subplot(212)
        plt.xticks(rotation=30)
        plt.pcolormesh(X,Y,los_synth_interp,cmap=cmap,vmin=zlims[0],vmax=zlims[1])
        plt.title('LOS modeled(m)')
        #plt.axis('equal')
        plt.colorbar()
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.grid()
        plt.suptitle(title)
        
    else:    
        ax=plt.subplot(211)
        ax.tick_params(labelbottom='off') 
        plt.scatter(lon,lat,c=los_data,cmap=cmap,vmin=zlims[0],vmax=zlims[1],s=50,lw=0)
        plt.title('LOS observed (m)')
        plt.colorbar()
        plt.grid()
        plt.ylabel('Latitude')
        #replt.axis('equal')
        plt.subplot(212)
        plt.xticks(rotation=30)
        plt.scatter(lon,lat,c=los_synth,cmap=cmap,vmin=zlims[0],vmax=zlims[1],s=50,lw=0)
        plt.title('LOS modeled(m)')
        #plt.axis('equal')
        plt.colorbar()
        plt.ylabel('Latitude')
        plt.xlabel('Longitude')
        plt.grid()
        plt.suptitle(title)

    if show:
        plt.show()
    
def tsunami_synthetics(home,project_name,run_name,run_number,gflist,t_lim,sort,scale):
    '''
    Plot synthetics vs real data
    
    gflist: The GF control fiel that decides what to plot/not plot
    datapath
    '''
    from obspy import read
    from numpy import genfromtxt,where,argsort,array,size
    import matplotlib.pyplot as plt
    import matplotlib
    from mudpy.green import stdecimate 
    from mudpy.forward import lowpass as lfilter
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='U')    
    lon=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[6],dtype='f')
    
    if size(sta)==1:
        sta=array([sta])
        lon=array([lon])
        lat=array([lat])
        gf=array([gf])
    
    
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    #Decide on sorting
    i=where(gf==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    nsta=len(i)
    
    fig, axarr = plt.subplots(nsta,1)  
    if len(i)==1:
        axarr=array([axarr])
    for k in range(len(i)):
        #tsun=read(datapath+sta[i[k]]+'.tsun')
        tsun=read(datapath+sta[i[k]]+'.sac')
        # tsun=read(datapath+'D43413.sac')
        tsun_synth=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.tsun')
        
        dt=tsun_synth[0].stats.starttime - tsun[0].stats.starttime
        
        #Make plot
        ax=axarr[k]
        ax.plot(tsun[0].times()/60-dt/60,tsun[0].data/scale,'k',tsun_synth[0].times()/60,tsun_synth[0].data/scale,'r')
        ax.grid(which='both')
        ax.yaxis.set_ticklabels([])
        ax.set_xlim(t_lim)
        ax.yaxis.grid(False)
        ax.yaxis.set_ticks([])
        #Annotations
        trange=t_lim[1]-t_lim[0]
        sign=1.
        if abs(min(tsun[0].data))>max(tsun[0].data):
            sign=-1. 
        tsun_max='%.3f' % (sign*max(abs(tsun[0].data))/scale)
        sign=1.
        if abs(min(tsun_synth[0].data))>max(tsun_synth[0].data):
            sign=-1. 
        tsun_synth_max='%.3f' % (sign*max(abs(tsun_synth[0].data))/scale)
        sign=1.
        tsun_lims=ax.get_ylim()
        tsun_range=tsun_lims[1]-tsun_lims[0]
        
        ax.annotate(tsun_max,xy=(t_lim[1]-0.2*trange,tsun_lims[0]+0.02*tsun_range),fontsize=12)
        ax.annotate(tsun_synth_max,xy=(t_lim[1]-0.2*trange,tsun_lims[0]+0.7*tsun_range),fontsize=12,color='red')
        #Station name
        ax.set_ylabel(sta[i[k]],rotation=90)
        #axn.set_title('North (m)')
        if k!=len(i)-1:
            ax.xaxis.set_ticklabels([])
        #    xtick=ax.xaxis.get_majorticklocs()
        #    ix=[1,3,5]
        #    xtick=xtick[ix]
        #    xticklabel=['','50','','150','','250','']
        if k==len(i)-1: #Last plot
            ax.set_xlabel('Minutes after Origin Time')
            ax.set_title('Tsunami (cm)')
            #ax.xaxis.set_ticklabels(xticklabel)
            #axn.xaxis.set_ticks(xtick)
            #axe.xaxis.set_ticks(xtick)
            #axu.xaxis.set_ticks(xtick)
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.8, top=0.85, wspace=0, hspace=0)
    plt.show()
                

def ABIC(home,project_name,run_name):
    '''
    plot values of ABIC vs smoothing parameter for model selection
    '''
    from glob import glob
    from numpy import zeros,argmin
    import matplotlib.pyplot as pl
    
    #Text rendering
    pl.rc('font',family='serif')
    #Get list of log files
    outdir=home+project_name+'/output/inverse_models/models/'
    plotdir=home+project_name+'/plots/'
    logs=glob(outdir+'*'+run_name+'.????.log')
    ABIC=zeros(len(logs))
    ls=zeros(len(logs))
    print('Gathering statistics for '+str(len(logs))+' inversions...')
    for k in range(len(logs)):
        with open(logs[k]) as f:
            for line in f:
                if 'ABIC' in line:
                    ABIC[k]=float(line.split('=')[1])
                if 'lambda_spatial' in line:
                    ls[k]=float(line.split('=')[1])
    #Get the minimum
    imin=argmin(ABIC)
    #Plot the thing
    pl.figure()
    pl.semilogx(ls,ABIC,'k',linewidth=2)
    pl.grid(which='both')
    pl.semilogx(ls[imin],ABIC[imin],marker='*',color='r',markersize=14)
    pl.xlabel(r'$\lambda$',fontsize=14)
    pl.ylabel('ABIC',fontsize=14)
    pl.annotate(r'$\lambda$'+'='+str(ls[imin]),xy=(ls[imin],ABIC[imin]),xytext=(ls[imin],ABIC[imin]-0.05*(max(ABIC)-ABIC[imin])))
    pl.title('Run Name: '+run_name)
    pl.savefig(plotdir+run_name+'.ABIC.png')
    print('ABIC is minimized at inversion '+logs[imin])
    print('... lambda = '+repr(ls[imin]))
    pl.show()
    
def ABIC2D(home,project_name,run_name,ABICmin,ABICmax):
    '''
    plot 2D values of ABIC vs smoothing parameter for model selection
    '''
    from glob import glob
    from numpy import zeros,argmin,log10,linspace
    import matplotlib.pyplot as plt
    from matplotlib import mlab as ml
    
    #Text rendering
    plt.rc('font',family='serif')
    #Get list of log files
    outdir=home+project_name+'/output/inverse_models/models/'
    if type(run_name)==list:
        for k in range(len(run_name)):
            if k==0:
                logs=glob(outdir+'*'+run_name[k]+'.????.log')
            else:
                logs+=glob(outdir+'*'+run_name[k]+'.????.log') 
    else:
        logs=glob(outdir+'*'+run_name+'.????.log')
    ABIC=zeros(len(logs))
    ls=zeros(len(logs))
    lt=zeros(len(logs))
    print('Gathering statistics for '+str(len(logs))+' inversions...')
    for k in range(len(logs)):
        with open(logs[k]) as f:
            for line in f:
                if 'ABIC' in line:
                    ABIC[k]=float(line.split('=')[1])
                if 'lambda_spatial' in line:
                    ls[k]=log10(float(line.split('=')[1]))
                if 'lambda_temporal' in line:
                    lt[k]=log10(float(line.split('=')[1]))
    #Get the minimum
    imin=argmin(ABIC)
    #Grid
    ABIC=ABIC/1000
    lsi=linspace(ls.min(),ls.max(),100)
    lti=linspace(lt.min(),lt.max(),100)
    ABICi=ml.griddata(ls,lt,ABIC,lsi,lti,interp='linear')
    #Plot the thing
    plt.figure()
    plt.pcolormesh(lsi,lti,ABICi,cmap=plt.cm.spectral_r,vmin=ABICmin,vmax=ABICmax)
    cb=plt.colorbar()
    plt.scatter(ls,lt,c='w',marker='o',s=30)
    plt.ylabel(r'$\log(\lambda_t$)',fontsize=18)
    plt.xlabel(r'$\log(\lambda_s$)',fontsize=18)
    cb.set_label(r'ABIC $(\times10^3)$')
    plt.scatter(ls[imin],lt[imin],marker='*',color='r',s=125)
    plt.title(r'$\lambda_s^{min} = $%.4e , $\lambda_t^{min} = $%.4e' % (10**ls[imin],10**lt[imin]))
    plt.show()
    print('ABIC is minimized at inversion '+logs[imin])
    print('... ls = '+repr(10**ls[imin])+' , lt = '+repr(10**lt[imin]))
        
    
def coherence(home,project_name,run_name,run_number,GF_list,vord,sort,f_lims):
    '''
    Plot coherences
    '''
    import matplotlib.pyplot as plt
    from numpy import load,where,genfromtxt,array,log10,argsort
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    lon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[2],dtype='f')
    datapath=home+project_name+'/analysis/frequency/'
    if vord.lower()=='d':
        kgf=0 #disp
        suffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        suffix='vel'
    i=where(gf[:,kgf]==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    #Initalize canvas
    fig, axarr = plt.subplots(len(i), 3)  
    for k in range(len(i)):
        #Read coherences
        coh=load(datapath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.coh.npz')
        fn=coh['fn']
        fe=coh['fe']
        fu=coh['fu']
        cn=coh['cn']
        ce=coh['ce']
        cu=coh['cu']
        #Let's plot them
        #get current axis
        axn=axarr[k,0]
        axe=axarr[k,1]
        axu=axarr[k,2]
        #Plot
        axn.semilogx(fn,cn)
        axe.semilogx(fe,ce,'g')
        axu.semilogx(fu,cu,'r')
        axn.grid(which='both')
        axe.grid(which='both')
        axu.grid(which='both')
        #Arrange axes
        axn.set_xlim(f_lims)
        axe.set_xlim(f_lims)
        axu.set_xlim(f_lims)
        axn.set_ylim([0,1])
        axe.set_ylim([0,1])
        axu.set_ylim([0,1])
        #Text labels
        axn.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axe.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axu.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axe.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        axn.yaxis.set_ticklabels(['','0.5','','1.0'])
        if k!=len(i)-1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
        if k==0: #First plot add some labels
            axn.set_title('North')
            if vord.lower()=='d':
                axe.set_title('Displacement Coherence\nEast')
            else:
                axe.set_title('Velocity Coherence\nEast')
            axu.set_title('Up')
        if k==len(i)-1: #Last plot
            axe.set_xlabel('Frequency (Hz)')
            #l=axe.get_xticks().tolist()
            #l[0]=''
            #axe.xaxis.set_ticklabels(l)
            #l=axu.get_xticks().tolist()
            #l[0]=''
            #axu.set_xticklabels(l)
        #Annotate with station name
        xyannot=(axn.get_xlim()[0]+0.01*log10((log10(axn.get_xlim()[1])-log10(axn.get_xlim()[0]))),axn.get_ylim()[0]+0.05)
        axn.annotate(sta[i[k]], xy=xyannot)
    plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, wspace=0, hspace=0)
    
def psds(home,project_name,run_name,run_number,GF_list,vord,sort,f_lims):
    '''
    Plot coherences
    '''
    import matplotlib.pyplot as plt
    from numpy import load,where,genfromtxt,array,log10,argsort
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    lon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1],dtype='f')
    lat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[2],dtype='f')
    datapath=home+project_name+'/analysis/frequency/'
    if vord.lower()=='d':
        kgf=0 #disp
        suffix='kdisp'
    elif vord.lower()=='v':
        kgf=1 #disp
        suffix='kvel'
    i=where(gf[:,kgf]==1)[0]  
    if sort.lower()=='lon':
        j=argsort(lon[i])[::-1]
        i=i[j]
    elif sort.lower()=='lat':
        j=argsort(lat[i])[::-1] 
        i=i[j]
    #Initalize canvas
    fig, axarr = plt.subplots(len(i), 3)  
    for k in range(len(i)):
        #Read coherences
        psd=load(datapath+sta[i[k]]+'.'+suffix+'.psd.npz')
        fn=psd['fn']
        fe=psd['fe']
        fu=psd['fu']
        npsd=psd['npsd']
        epsd=psd['epsd']
        upsd=psd['upsd']
        #Let's plot them
        #get current axis
        axn=axarr[k,0]
        axe=axarr[k,1]
        axu=axarr[k,2]
        #Plot
        axn.semilogx(fn,npsd)
        axe.semilogx(fe,epsd,'g')
        axu.semilogx(fu,upsd,'r')
        axn.grid(which='both')
        axe.grid(which='both')
        axu.grid(which='both')
        #Arrange axes
        axn.set_xlim(f_lims)
        axe.set_xlim(f_lims)
        axu.set_xlim(f_lims)
        axn.set_ylim([npsd.min(),npsd.max()])
        axe.set_ylim([npsd.min(),npsd.max()])
        axu.set_ylim([npsd.min(),npsd.max()])
        #Text labels
        #axn.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        #axe.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        #axu.yaxis.set_ticks(array([0.25,0.5,0.75,1.0]))
        axe.yaxis.set_ticklabels([])
        axn.yaxis.set_ticklabels([])
        axu.yaxis.set_ticklabels([])
        #axn.yaxis.set_ticklabels(['','0.5','','1.0'])
        if k!=len(i)-1:
            axn.xaxis.set_ticklabels([])
            axe.xaxis.set_ticklabels([])
            axu.xaxis.set_ticklabels([])
        if k==0: #First plot add some labels
            axn.set_title('North')
            if vord.lower()=='d':
                axe.set_title('Displacement PSD(dB)\nEast')
            else:
                axe.set_title('Velocity PSD(dB)\nEast')
            axu.set_title('Up')
        if k==len(i)-1: #Last plot
            axe.set_xlabel('Frequency (Hz)')
            #l=axe.get_xticks().tolist()
            #l[0]=''
            #axe.xaxis.set_ticklabels(l)
            #l=axu.get_xticks().tolist()
            #l[0]=''
            #axu.set_xticklabels(l)
        #Annotate with station name
        xyannot=(axn.get_xlim()[0]+0.01*log10((log10(axn.get_xlim()[1])-log10(axn.get_xlim()[0]))),axn.get_ylim()[0]+0.05)
        axn.annotate(sta[i[k]], xy=xyannot)
    plt.subplots_adjust(left=0.25, bottom=0.05, right=0.75, top=0.9, wspace=0, hspace=0)
            
def average_coherence(home,project_name,run_name,run_number,GF_list,vord,num_components):
    '''
    Plot coherences
    '''
    import matplotlib.pyplot as plt
    from numpy import load,where,genfromtxt,zeros,interp
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/analysis/frequency/'
    if vord.lower()=='d':
        kgf=0 #disp
        suffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        suffix='vel'
    i=where(gf[:,kgf]==1)[0]  
    for k in range(len(i)):
        #Read coherences
        coh=load(datapath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.coh.npz')
        fn=coh['fn']
        fe=coh['fe']
        fu=coh['fu']
        cn=coh['cn']
        ce=coh['ce']
        cu=coh['cu'] 
        if k==0:
            f=fn
            c1=zeros(cn.shape)
            c1e=zeros(cn.shape)
            c1n=zeros(cn.shape) 
            c1u=zeros(cn.shape)    
        if num_components==1: #Average all
            try:
                c1+=cn
                c1+=ce
                c1+=cu
            except: #Coherence is the wrong size
                cn=interp(f,fn,cn)
                ce=interp(f,fn,ce)
                cu=interp(f,fn,cu)
                c1+=cn
                c1+=ce
                c1+=cu
        else:
            try:
                c1e+=ce
                c1n+=cn
                c1u+=cu
            except: #Coherence is the wrong size
                cn=interp(f,fn,cn)
                ce=interp(f,fn,ce)
                cu=interp(f,fn,cu)
                c1n+=cn
                c1e+=ce
                c1u+=cu
    #Normalize
    c1=c1/(3*len(i))
    c1e=c1e/len(i)
    c1n=c1n/len(i)
    c1u=c1u/len(i)
    #Plot
    plt.figure()
    if num_components==1:
        plt.semilogx(fn,c1)
        plt.fill_between(fn,y1=0,y2=c1,color='r',alpha=0.5)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Mean Coherence')
        plt.grid(which='both')
        plt.xlim(f.min(),f.max())
    else:
        plt.semilogx(fn,c1n,fe,c1e,fu,c1u)
        plt.legend(['North','East','Up'])
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Mean Coherence')
        plt.grid(which='both')
        plt.xlim(f.min(),f.max())
        
def stf_spectra(home,project_name,rupt,flims,dlims,normalize=True,stacks=None):
    '''
    Plot source time functions
    '''
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from numpy import load,genfromtxt,unique,arange,where,intersect1d
    
    
    datapath=home+project_name+'/analysis/frequency/'
    f=genfromtxt(rupt)
    rupt=rupt.split('/')[-1]
    depth=f[:,3]
    num=f[:,0]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    nfault=len(unum)
    depth=depth[0:nfault]
    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0,0],[0,0]]
    cm = plt.get_cmap('brg') 
    levels = arange(dlims[0],dlims[1]+0.1,0.1)
    plt.figure()
    c = plt.contourf(Z, levels, cmap=cm)
    plt.clf()
    cNorm  = colors.Normalize(vmin=depth.min(), vmax=depth.max())
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    if stacks==None:
        for k in range(nfault):
            #Read coherences
            P=load(datapath+rupt.split('.')[0]+'.'+rupt.split('.')[1]+'.sub'+str(k).rjust(4,'0')+'.stfpsd.npz')
            f=P['freq']
            psd=P['psd']
            i1=where(f>=flims[0])[0]
            i2=where(f<=flims[1])[0]
            i=intersect1d(i1,i2)
            f=f[i]
            psd=psd[i]
            if normalize==True:
                psd=psd/psd.max()
            colorVal = scalarMap.to_rgba(depth[k])
            plt.loglog(f,psd,color=colorVal)
    else:
        for k in range(len(stacks)):
            current=stacks[k]
            for ks in range(len(current)):
                #Read coherences
                P=load(datapath+rupt.split('.')[0]+'.'+rupt.split('.')[1]+'.sub'+str(current[ks]).rjust(4,'0')+'.stfpsd.npz')
                f=P['freq']
                psd=P['psd']
                i1=where(f>=flims[0])[0]
                i2=where(f<=flims[1])[0]
                i=intersect1d(i1,i2)
                f=f[i]
                psd=psd[i]
                if ks==0:
                    p=psd
                    z=depth[current[ks]]
                else:
                    p+=psd
                    z+=depth[current[ks]]
            if normalize==True:
                p=p/p.max()
            else:
                p=p/len(current)
            z=z/len(current)
            colorVal = scalarMap.to_rgba(z)
            plt.loglog(f,p,color=colorVal)
    plt.xlim(flims)
    plt.grid(which='both')
    cb=plt.colorbar(c)
    cb.set_label('Subfault depth (km)')
    plt.xlabel('Frequency (Hz)')
    if normalize==True:
        plt.ylabel('Normalized Power')
    else:
        plt.ylabel('PSD (m/s)'+r'$^2$'+'/Hz')
    plt.title('Slip Rate Spectra')


def dtopo_slices(dtopo_file,fault_file,out):
    '''
    Plot dtopo file
    '''
    from numpy import genfromtxt,unique,where,meshgrid
    from scipy.interpolate import griddata
    from matplotlib import pyplot as plt
    
    dtopo=genfromtxt(dtopo_file)
    f=genfromtxt(fault_file)
    t=unique(dtopo[:,0])
    #Get z ranges
    z1=dtopo[:,3].min()
    z2=dtopo[:,3].max()
    zrange=max([-z1,z2])
    #Loop over time slices
    for kt in range(len(t)):
        print('t = '+str(t[kt]))
        i=where(dtopo[:,0]==t[kt])[0]
        lon=dtopo[i,1]
        lat=dtopo[i,2]     
        Lon=unique(dtopo[i,1])
        Lat=unique(dtopo[i,2])
        Lon,Lat=meshgrid(Lon,Lat)
        z=dtopo[i,3]
        Z=griddata((lon,lat),z,(Lon,Lat),method='linear')
        plt.figure();
        plt.imshow(Z,origin='lower',extent=(lon.min(),lon.max(),lat.min(),lat.max()),
                    vmin=-zrange,vmax=zrange,cmap=plt.cm.seismic);
        cb=plt.colorbar()
        cb.set_label('Vertical deformation (m)')
        plt.title('t = '+str(t[kt]))
        plt.xlabel('Longitude (deg)')
        plt.ylabel('Latitude (deg)')
        plt.scatter(f[:,1],f[:,2],c='k',marker='x')
        plt.savefig(out+str(kt).rjust(4,'0')+'vert.png')
        plt.close("all")
        

def plot_grd(grdfile,zlims,cmap,flip_lon=False,return_data=False,z_variable=None,contours=None,title=None):
    '''
    Quick plot of any GMT grd file, this will only work for NETCDF4 files, 
    i.e. if you use GMT5. If you are outdated and use NETCDF3 you can edit this
    to use scipy.io.netcdf_file instead.
    
    grdfile - path to file
    '''
    from netCDF4 import Dataset
    from numpy import meshgrid,genfromtxt
    import matplotlib.pyplot as plt

    grd = Dataset(grdfile, 'r', format='NETCDF4')
    if z_variable==None:
        try:
            x=grd.variables['x'][:]
            y=grd.variables['y'][:]
            z=grd.variables['z'][:]
        except:
            x=grd.variables['lon'][:]
            y=grd.variables['lat'][:]
            z=grd.variables['z'][:]
    else:
            x=grd.variables['lon'][:]
            y=grd.variables['lat'][:]
            z=grd.variables[z_variable][:]
    if flip_lon:
        x=x-360
    X,Y=meshgrid(x,y)
    plt.figure()
    if title==None:
        plt.title(grdfile)
    else:
        plt.title(title)
    plt.pcolormesh(X,Y,z,vmin=zlims[0],vmax=zlims[1],cmap=cmap)
    plt.colorbar()
    
    if contours != None:
        plt.contour(X,Y,z,levels=contours)
    
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()
    
    if return_data:
        return X,Y,z




def plot_velmod(vel_mod_file1,vel_mod_file2=None,zmax=60,label1='Cascadia',label2='GIL7'):
    '''
    Plot a velocity model
    '''
    
    from numpy import genfromtxt,zeros,arange
    from matplotlib import pyplot as plt
    
    v1=genfromtxt(vel_mod_file1)
    #Init
    pvel1=zeros(len(v1)*2)
    svel1=zeros(len(v1)*2)
    rho1=zeros(len(v1)*2)
    Qp1=zeros(len(v1)*2)
    Qs1=zeros(len(v1)*2)
    z1=zeros(len(v1)*2)
    #interleave
    i1=arange(0,len(z1),2)
    i2=arange(1,len(z1),2)
    pvel1[i1]=v1[:,2] ; pvel1[i2]=v1[:,2]
    svel1[i1]=v1[:,1] ; svel1[i2]=v1[:,1]
    rho1[i1]=v1[:,3] ; rho1[i2]=v1[:,3]
    Qp1[i1]=v1[:,5] ; Qp1[i2]=v1[:,5]
    Qs1[i1]=v1[:,4] ; Qs1[i2]=v1[:,4]
    z1[0]=0
    z1[1]=v1[0,0]
    for k in range(1,len(v1)-1):
        z1[2*k]=z1[2*k-1]
        z1[2*k+1]=z1[2*k]+v1[k,0]
    z1[-2]=z1[-3]
    z1[-1]=zmax
    mu1=rho1*1000*(svel1*1000)**2
    
    if vel_mod_file2!=None:
        v2=genfromtxt(vel_mod_file2)
        #Init
        pvel2=zeros(len(v2)*2)
        svel2=zeros(len(v2)*2)
        rho2=zeros(len(v2)*2)
        Qp2=zeros(len(v2)*2)
        Qs2=zeros(len(v2)*2)
        z2=zeros(len(v2)*2)
        #interleave
        i1=arange(0,len(z2),2)
        i2=arange(1,len(z2),2)
        pvel2[i1]=v2[:,2] ; pvel2[i2]=v2[:,2]
        svel2[i1]=v2[:,1] ; svel2[i2]=v2[:,1]
        rho2[i1]=v2[:,3] ; rho2[i2]=v2[:,3]
        Qp2[i1]=v2[:,5] ; Qp2[i2]=v2[:,5]
        Qs2[i1]=v2[:,4] ; Qs2[i2]=v2[:,4]
        z2[0]=0
        z2[1]=v2[0,0]
        for k in range(1,len(v2)-1):
            z2[2*k]=z2[2*k-1]
            z2[2*k+1]=z2[2*k]+v2[k,0]
        z2[-2]=z2[-3]
        z2[-1]=zmax
        mu2=rho2*1000*(svel2*1000)**2
        
    #plotaroo
    plt.figure(figsize=(20,5))
    
    
    plt.subplot(161)
    plt.plot(pvel1,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(pvel2,-z2,'r',lw=2)
        plt.legend([label1,label2],loc=3)
    plt.xticks(rotation=-45)
    plt.xlabel('Vp (km/s)')
    plt.ylabel('Depth (km)')
    plt.show()
        
    plt.subplot(162)
    plt.plot(svel1,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(svel2,-z2,'r',lw=2)
    plt.xticks(rotation=-45)
    plt.tick_params(axis='y', labelleft='off')
    plt.xlabel('Vs (km/s)')
    
    plt.subplot(163)
    plt.plot(rho1,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(rho2,-z2,'r',lw=2)
    plt.xticks(rotation=-45)
    plt.tick_params(axis='y', labelleft='off')
    plt.xlabel('Density (g/cm^3)')
    
    plt.subplot(164)
    plt.plot(Qp1,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(Qp2,-z2,'r',lw=2)
    plt.xticks(rotation=-45)
    plt.tick_params(axis='y', labelleft='off')
    plt.xlabel('Qp')
    
    plt.subplot(165)
    plt.plot(Qs1,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(Qs2,-z2,'r',lw=2)
    plt.xticks(rotation=-45)
    plt.tick_params(axis='y', labelleft='off')
    plt.xlabel('Qs') 
    
    plt.subplot(166)
    plt.plot(mu1/1e9,-z1,'k',lw=2)
    if vel_mod_file2!=None:
        plt.plot(mu2/1e9,-z2,'r',lw=2)
    plt.xticks(rotation=-45)
    plt.tick_params(axis='y', labelleft='off')
    plt.xlabel('Rigidity (GPa)')    
    
    plt.subplots_adjust(top=0.95,bottom=0.2,right=0.95)
    plt.show()

    

def plot_mesh(mesh,lw=0.5,zaspect=0.2):
    
    """
    Plots a mesh given as a NumPy array of vertices.
    
    Args:
        mesh (numpy.ndarray): A NumPy array of vertices, where each row contains
            the coordinates of three vertices that form a triangle.
        lw (float, optional): The line width of the mesh. Default is 0.5.
        zaspect (float, optional): The aspect ratio of the z-axis relative to the
            x- and y-axes. Default is 0.2.
    
    Returns:
        None
    
    Example:
        >>> import numpy as np
        >>> mesh = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
        >>> plot_mesh(mesh)
    """
        
    
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    for k in range(len(mesh)):
        
        x = [mesh[k,4],mesh[k,7],mesh[k,10],mesh[k,4]]
        y = [mesh[k,5],mesh[k,8],mesh[k,11],mesh[k,5]]
        z = [mesh[k,6],mesh[k,9],mesh[k,12],mesh[k,6]]
        
        plt.plot(x,y,z,lw=lw,c='b')
        
    ax.set_box_aspect((1,1,zaspect))

    ax.set_xlabel('lon')
    ax.set_ylabel('lat')
    ax.set_zlabel('depth (km)')





#########                  Supporting tools                       ##############

def slip2geo(ss,ds,strike):
    '''
    Determine geogrpahical orientation of rake vector
    '''
    from numpy import deg2rad,sin,cos
    
    #Normalize slips
    ds=ds/((ds**2+ss**2)**0.5)
    ss=ss/((ds**2+ss**2)**0.5)
    #determine contribution of ds and ss slips
    xds=ds*sin(deg2rad(strike-90))
    yds=ds*cos(deg2rad(strike-90))
    xss=ss*sin(deg2rad(strike))
    yss=ss*cos(deg2rad(strike))
    #Add em up
    x=xss+xds
    y=yss+yds
    return x,y
