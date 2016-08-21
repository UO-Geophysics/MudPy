'''
D.Melgar
04/2014

Some routines to make quick and not so quick plots of the forward modeling and 
inversion results
'''

import matplotlib
from matplotlib import cm

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
whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

'''
numpy.interp(x, xp, fp, left=None, right=None)[source]
One-dimensional linear interpolation.

Returns the one-dimensional piecewise linear interpolant to a function with given values at discrete data-points.

Parameters:	
    x : array_like
    The x-coordinates of the interpolated values.
    
    xp : 1-D sequence of floats
    The x-coordinates of the data points, must be increasing.
    
    fp : 1-D sequence of floats
    The y-coordinates of the data points, same length as xp.
    
    left : float, optional
    Value to return for x < xp[0], default is fp[0].

    right : float, optional
    Value to return for x > xp[-1], default is fp[-1].
    
Returns:	
    y : {float, ndarray}
    The interpolated values, same shape as x.

Raises:	

    ValueError :
    If xp and fp have different length
'''


def quick_model(rupt):
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
    #Plot
    plt.figure()
    plt.scatter(lon,lat,marker='o',c=slip,s=250,cmap=whitejet)
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.quiver(lon,lat,x,y,color='green',width=0.0013)
    plt.grid()
    plt.title(rupt)
    plt.show()
    
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


def slip3D(rupt,marker_size=60,clims=None):
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

    #Plot it
    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(111, projection='3d')
    if clims==None:
        p=ax.scatter(lon, lat, depth, c=slip,cmap=whitejet, marker='o',s=marker_size,lw=0)
    else:
        p=ax.scatter(lon, lat, depth, c=slip,cmap=whitejet, marker='o',s=marker_size,vmin=clims[0],vmax=clims[1],lw=0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_zlabel('Depth (km)')
    cb=fig.colorbar(p)
    cb.set_label('Slip (m)')
    plt.subplots_adjust(left=0.1, bottom=0.1, right=1.0, top=0.9, wspace=0, hspace=0)
    plt.title(rupt)
    plt.show()


def plot_insar(home,project_name,GF_list,(los_min,los_max)):
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


def tile_slip(rupt,nstrike,ndip,(slip_bounds),geographic=False,epicenter=0,epicenter_line=0,thresh=0):
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
    from obspy.core.util.geodetics import gps2DistAzimuth
    
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
    lon_afters=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=4)
    lat_afters=genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=3)
    depth_afters=-genfromtxt('/Users/dmelgar/Lefkada2015/afters/aftershocks_NOA_reloc.txt',usecols=5)
    if geographic==True: #Get geographic coordinates to compute along strike and along dip distance
        lon=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,1] #Only compute line at the epicenter depth
        lat=f[(epicenter_line-1)*nstrike:epicenter_line*nstrike,2]    
        depth=-f[:,3]
        depth=depth[0:len(unum)]
        along_strike=zeros(nstrike)
        for k in range(len(lat)):
            out=gps2DistAzimuth(epicenter[1],epicenter[0],lat[k],lon[k])
            if lat[k]<epicenter[1]: #It's to the south
                along_strike[k]=-out[0]/1000
            else:
                along_strike[k]=out[0]/1000
        #Now tile
        along_strike=tile(along_strike,ndip)
        #Process the aftershocks
        along_strike_afters=zeros(len(lon_afters))
        for k in range(len(lat_afters)):
            out=gps2DistAzimuth(epicenter[1],epicenter[0],lat_afters[k],lon_afters[k])
            if lat_afters[k]<epicenter[1]: #It's to the south
                along_strike_afters[k]=-out[0]/1000
            else:
                along_strike_afters[k]=out[0]/1000
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
        plt.scatter(istrike,idip,marker='o',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
        plt.ylabel('Along-dip index')
        plt.xlabel('Along-strike index')
        plt.xlim(istrike.min()-1,istrike.max()+1)
        plt.ylim(idip.min()-1,idip.max()+1)
        plt.quiver(istrike,idip,rakess,rakeds,color='green',width=0.0001)
        plt.axis('equal')
        plt.grid()    
        plt.title(rupt)
    else:
        rakess=rakess*slip
        rakeds=rakeds*slip
        plt.figure(num=None, figsize=(13, 3.5), dpi=80)
        #plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        plt.scatter(along_strike,depth,marker='s',linewidth=0.5,edgecolor='#CCCCCC',c=slip,s=250,cmap=plt.cm.magma_r,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
        #plt.scatter(along_strike_afters,depth_afters,marker='.',c='#404040',s=35)
        plt.ylabel('Depth (km)')
        plt.xlabel('Along-strike distance (km)')
        #plt.xlim(along_strike.min()-5,along_strike.max()+5)
        #plt.ylim(depth.min()-5,depth.max()+5)   
        plt.xlim(-39,41)
        plt.ylim(-17,1.5)       
        plt.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='#00FF00',s=350,linewidth=2)
        for k in range(len(along_strike)):
            scale_slip=slip[k]/slip.max()
            plt.quiver(along_strike[k],depth[k],rakess[k]/sqrt(rakess[k]**2+rakeds[k]**2),rakeds[k]/sqrt(rakess[k]**2+rakeds[k]**2),color='green',width=0.002,scale=50/scale_slip)
    plt.annotate('North',xy=(28,0),fontsize=16)
    plt.annotate('South',xy=(-36,0),fontsize=16)
    #plt.title(r'2015 Lefkada $M_w6.55$, $v_r=2.6$km/s, $\sigma=020$, $\delta=65$',fontsize=16)
    cb.set_label('Slip(m)')
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.92, top=0.92, wspace=0, hspace=0)
    plt.show()


def tile_resolution(rupt,resfile,nstrike,ndip,(res_min,res_max),epicenter=0,epicenter_line=0):
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


   
def tile_slip_movieframes(home,project_name,sliprate_path,nstrike,ndip,(slip_min,slip_max),dt,geographic=False,epicenter=0,epicenter_line=0):
    '''
    Quick and dirty plot of a .rupt file
    epicenter is the coordinates, epcienter line is the down dip lien number where 
    the epcienter is
    '''
    
    from numpy import genfromtxt,zeros,tile,linspace,pi,cos,sin,ones,meshgrid,arange,r_
    import matplotlib.pyplot as plt
    from obspy.core.util.geodetics import gps2DistAzimuth
    from glob import glob
    from string import rjust,ljust
    from matplotlib.mlab import griddata
    
    #Get sliprate files
    files=glob(sliprate_path+'*.sliprate')
    for kframe in range(len(files)):
        print kframe
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
        t=dt*kframe #Current time
        print "t="+str(t)
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
        x=linspace(along_strike.min()-0.5,along_strike.max()+0.5,100)  #Adjsut for the width of the subfault
        y=linspace(depth.min()-0.49,depth.max()+0.49,100)   #Adjsut for height of subfault
        X, Y = meshgrid(x, y)
        sliprate=griddata(along_strike,depth,slip,X,Y,interp='linear')
        plt.scatter(along_strike,depth,marker='s',linewidth=0.0,edgecolor='',c=slip,s=250,cmap=whitejet,vmin=slip_min,vmax=slip_max)
        cb=plt.colorbar()
        plt.contourf(X,Y,sliprate,100,vmin=slip_min,vmax=slip_max,cmap=whitejet)
        plt.grid()
        #End interpolated 
        plt.ylabel('Depth (km)')
        plt.xlabel('Along-strike distance (km)')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.scatter(0,-epicenter[2],marker='*',edgecolor='k',facecolor='#00FF00',s=350,linewidth=2)
        #plt.scatter(xcent,zcent,marker='D',edgecolor='black',facecolor='',s=120,linewidth=2)

        plt.scatter(xaf,zaf,edgecolor='k',s=5)
        plt.plot(x15,y15,'--',c='grey')
        plt.plot(x20,y20,'--',c='grey')
        plt.plot(x25,y25,'--',c='grey')

        cb.set_label('Slip rate (m/s)')
        plot_name=home+project_name+'/plots/sliprate.'+rjust(str(kframe),4,'0')+'.png'
        plt.title('t = '+ljust(str(kframe*dt),4,'0')+'s')
        plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.9, wspace=0, hspace=0)
        plt.savefig(plot_name)     
        ts=kframe*dt
        if ts==1.0 or ts==2.0 or ts==3.0 or ts==4.0 or ts==5.0 or ts==6.0 or ts==7.0 or ts==8.0:
            plot_name=home+project_name+'/plots/sliprate.'+rjust(str(kframe),4,'0')+'.pdf'
            plt.savefig(plot_name) 
            print 'Saved PDF frame'
        plt.close("all")
    
        
            
def panel_tile_slip(home,project_name,sliprate_path,nstrike,ndip,(slip_min,slip_max),nframes,geographic=False,epicenter=0,epicenter_line=0):
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
    from string import rjust,ljust
    from matplotlib.mlab import griddata
    from matplotlib import colorbar
    
    #Get sliprate files
    files=[]
    for k in range(len(nframes)):
        frame=rjust(str(nframes[k]),4,'0')
        files.append(sliprate_path+frame+'.slip')
    #Now intialize figure
    fig, axarr = plt.subplots(5, 2,figsize=(8,4))
    axarr=axarr.transpose()
    axarr=reshape(axarr,(10,))
    #loop through frames
    for kframe in range(len(files)):
        ax=axarr[kframe] #current axis
        print kframe
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
        print "t="+str(t)
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
                    
                            
                       
def tile_moment(rupt,epicenter,nstrike,ndip,covfile,beta=0,(vfast,vslow)=(0,0),shade=False):
    '''
    Tile plot of subfault source-time functions
    '''
    import matplotlib.pyplot as plt
    from matplotlib import cm
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
    iss=2*arange(0,len(all)/2,1)
    ids=2*arange(0,len(all)/2,1)+1
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
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
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
        t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        if covfile !=None:
            t1plus,M1plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_plus[0])
            t1minus,M1minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_minus[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            if covfile !=None:
                t2plus,M2plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_plus[kwin+1])
                t2minus,M2minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_minus[kwin+1])
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
        ax=axarr[int(idip[kfault]), int(istrike[kfault])]
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
        ax.plot(t1, M1,color='k')
        
        ax.grid()
        ax.set_xlim([t1[0],t1[-1]])
        ax.xaxis.set_ticks(linspace(t1[0],t1[-1],5))
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
    #Go back and rescale all subplots by maximum moment
    for k in range(ndip):
        for k2 in range(nstrike):
            ax=axarr[k,k2]
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
    print 'Maximum moment was '+str(Mmax)+'N-m'
    return tout,Mout


def source_time_function(rupt,epicenter):
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
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        if kfault==0:#Get first source time function
            t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
    #Get power
    exp=floor(log10(M1.max()))
    M1=M1/(10**exp)
    plt.figure()
    plt.fill(t1,M1,'b',alpha=0.5)
    plt.plot(t1,M1,color='k')
    plt.grid()
    plt.xlabel('Time(s)')
    plt.ylabel('Moment Rate ('+r'$\times 10^{'+str(int(exp))+r'}$Nm/s)')
    plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.7, wspace=0, hspace=0)
    return t1,M1
    
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
    from string import rjust
    
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
        print str(kslice)+'/'+str(len(tslice)-1)
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
            plt.savefig(out+rjust(str(kslice),4,'0')+'.kin_slice.png')
        else:
            cb.set_label('Cumulative Slip (m)')
            plt.title('t = '+str(tslice[kslice+1])+'s')
            plt.savefig(out+rjust(str(kslice),4,'0')+'.kin_cumulative.png')
        plt.close("all")
    

def plot_data(home,project_name,gflist,vord,decimate,lowpass,t_lim,sort,scale,k_or_g):
    '''
    Plot synthetics vs real data
    
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
      
def synthetics(home,project_name,run_name,run_number,gflist,vord,decimate,lowpass,t_lim,sort,scale,k_or_g,uncert=False):
    '''
    Plot synthetics vs real data
    
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
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    if vord.lower()=='d':
        kgf=0 #disp
        if k_or_g.lower()=='kal':
            datasuffix='kdisp'
        else:
            datasuffix='disp'
        synthsuffix='disp'
    elif vord.lower()=='v':
        kgf=1 #disp
        datasuffix='vel'
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
        ns=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.n.sac')
        es=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.e.sac')
        us=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.u.sac')
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
        axn.plot(n[0].times(),n[0].data,'k',ns[0].times(),ns[0].data,'r')
        axn.grid(which='both')
        axe.plot(e[0].times(),e[0].data,'k',es[0].times(),es[0].data,'r')
        axe.grid(which='both')
        axu.plot(u[0].times(),u[0].data,'k',us[0].times(),us[0].data,'r')
        if uncert==True:
            axn.fill_between(n[0].times(),n[0].data-0.01,n[0].data+0.01,alpha=0.2,color='k')
            axe.fill_between(e[0].times(),e[0].data-0.01,e[0].data+0.01,alpha=0.2,color='k')
            axu.fill_between(u[0].times(),u[0].data-0.03,u[0].data+0.03,alpha=0.2,color='k')
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
        if abs(min(ns[0].data))>max(ns[0].data):
            sign=-1. 
        nsmax='%.3f' % (sign*max(abs(ns[0].data)))
        sign=1.
        nlims=axn.get_ylim()
        nrange=nlims[1]-nlims[0]
        
        if abs(min(e[0].data))>max(e[0].data):
            sign=-1.         
        emax='%.3f' % (sign*max(abs(e[0].data)))
        sign=1.
        if abs(min(es[0].data))>max(es[0].data):
            sign=-1. 
        esmax='%.3f' % (sign*max(abs(es[0].data)))
        sign=1.
        elims=axe.get_ylim()
        erange=elims[1]-elims[0]
        
        if abs(min(u[0].data))>max(u[0].data):
            sign=-1. 
        umax='%.3f' % (sign*max(abs(u[0].data)))
        sign=1.
        if abs(min(us[0].data))>max(us[0].data):
            sign=-1 
        usmax='%.3f' % (sign*max(abs(us[0].data)))
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
        axn.annotate(nsmax,xy=(t_lim[1]-0.25*trange,nlims[0]+0.7*nrange),fontsize=12,color='red')
        axe.annotate(esmax,xy=(t_lim[1]-0.25*trange,elims[0]+0.7*erange),fontsize=12,color='red')
        axu.annotate(usmax,xy=(t_lim[1]-0.25*trange,ulims[0]+0.7*urange),fontsize=12,color='red')
        #Station name
        axn.set_ylabel(sta[i[k]],rotation=90)
        if k==0:
            if vord.lower()=='d':
                axn.set_title('North (m)')
                axe.set_title('East (m)')
                axu.set_title('Up (m)')
            else:
                axn.set_title('North (m/s)')
                axe.set_title('East (m/s)')
                axu.set_title('Up (m/s)')
        if k!=len(i)-1:
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
            #xticklabel=['','100','','200','','300'] #Maule preferred
            #xticklabel=['','','40','','80','','120','','160',''] #Iquique preferred
            #xticklabel=['','10','','30','','50',''] #Nepal preferred
            #xticklabel=['','5','','15','','25','','35',''] #Lefkada preferred
        if k==len(i)-1 and nsta>1: #Last plot
            axe.set_xlabel('Time (s)')
            #axn.xaxis.set_ticklabels(xticklabel)
            #axe.xaxis.set_ticklabels(xticklabel)
            #axu.xaxis.set_ticklabels(xticklabel)
            #axn.xaxis.set_ticks(xtick)
            #axe.xaxis.set_ticks(xtick)
            #axu.xaxis.set_ticks(xtick)
    #plt.subplots_adjust(left=0.2, bottom=0.05, right=0.8, top=0.95, wspace=0, hspace=0)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.8, top=0.85, wspace=0, hspace=0)

def static_synthetics(home,project_name,run_name,run_number,gflist,qscale):
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
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='S')
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
        neu=genfromtxt(datapath+sta[i[k]]+'.neu')
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
    for k in range(len(i)):
        plt.annotate(sta[i[k]],xy=(lon[k],lat[k]))
    plt.xlim([lon.min()-lonadjust,lon.max()+lonadjust])
    plt.ylim([lat.min()-latadjust,lat.max()+latadjust])
    plt.subplot(122)
    plt.quiver(lon,lat,zeros(len(u)),u,color='k',scale=qscale)
    plt.quiver(lon,lat,zeros(len(us)),us,color='r',scale=qscale)
    plt.grid()
    plt.title('Verticals')
    for k in range(len(i)):
        plt.annotate(sta[i[k]],xy=(lon[k],lat[k]))
    plt.xlim([lon.min()-lonadjust,lon.max()+lonadjust])
    plt.ylim([lat.min()-latadjust,lat.max()+latadjust])
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
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='S')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[7],dtype='f')
    datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    los_data=zeros(len(i))
    los_synth=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(datapath+sta[i[k]]+'.los')
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        los_data[k]=neu[0]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.los')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        los_synth[k]=neus[0]
    #Make plot
    out=c_[lon,lat,los_data,los_synth,los_data-los_synth]
    savetxt(home+project_name+'/analysis/'+run_name+'.'+run_number+'.insar.res',out,fmt='%.6f\t%.6f\t%8.5f\t%8.5f\t%8.5f',header='lon,lat,los_data(m),los_synthetic(m),data-synthetic(m)')
    plt.figure()
    plt.scatter(lon,lat,c=los_data-los_synth,cmap=matplotlib.cm.seismic,vmin=zlims[0],vmax=zlims[1],s=50)
    plt.title('LOS data - LOS predicted (m)')
    plt.colorbar()
    plt.grid()
    
    
def insar_results(home,project_name,run_name,run_number,gflist,zlims):
    '''
    Plot insar observed in one panel and insar modeled in the other
    
    gflist: The GF control file that decides what to plot/not plot
    datapath
    '''
    from numpy import genfromtxt,where,zeros,sqrt,c_,savetxt
    import matplotlib.pyplot as plt
    import matplotlib
    
    matplotlib.rcParams.update({'font.size': 14})
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=0,dtype='S')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[7],dtype='f')
    datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    los_data=zeros(len(i))
    los_synth=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(datapath+sta[i[k]]+'.los')
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        los_data[k]=neu[0]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.los')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        los_synth[k]=neus[0]
    #Make plot
    out=c_[lon,lat,los_data,los_synth,los_data-los_synth]
    savetxt(home+project_name+'/analysis/'+run_name+'.'+run_number+'.insar.res',out,fmt='%.6f\t%.6f\t%8.5f\t%8.5f\t%8.5f',header='lon,lat,los_data(m),los_synthetic(m),data-synthetic(m)')
    plt.figure()
    ax=plt.subplot(211)
    ax.tick_params(labelbottom='off') 
    plt.scatter(lon,lat,c=los_data,cmap=cm.jet,vmin=zlims[0],vmax=zlims[1],s=50,lw=0)
    plt.title('LOS observed (m)')
    plt.colorbar()
    plt.grid()
    plt.ylabel('Latitude')
    #replt.axis('equal')
    plt.subplot(212)
    plt.xticks(rotation=30)
    plt.scatter(lon,lat,c=los_synth,cmap=cm.jet,vmin=zlims[0],vmax=zlims[1],s=50,lw=0)
    plt.title('LOS modeled(m)')
    #plt.axis('equal')
    plt.colorbar()
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.grid()
            
def tsunami_synthetics(home,project_name,run_name,run_number,gflist,t_lim,sort,scale):
    '''
    Plot synthetics vs real data
    
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
    gf=genfromtxt(home+project_name+'/data/station_info/'+gflist,usecols=[6],dtype='f')
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
    for k in range(len(i)):
        tsun=read(datapath+sta[i[k]]+'.tsun')
        tsun_synth=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.tsun')
        #Make plot
        ax=axarr[k]
        ax.plot(tsun[0].times()/60,tsun[0].data/scale,'k',tsun_synth[0].times()/60,tsun_synth[0].data/scale,'r')
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
            #ax.xaxis.set_ticklabels(xticklabel)
            #axn.xaxis.set_ticks(xtick)
            #axe.xaxis.set_ticks(xtick)
            #axu.xaxis.set_ticks(xtick)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.8, top=0.85, wspace=0, hspace=0)
                

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
    print 'Gathering statistics for '+str(len(logs))+' inversions...'
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
    print 'ABIC is minimized at inversion '+logs[imin]
    print '... lambda = '+repr(ls[imin])
    pl.show()
    
def ABIC2D(home,project_name,run_name,(ABICmin,ABICmax)):
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
    print 'Gathering statistics for '+str(len(logs))+' inversions...'
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
    print 'ABIC is minimized at inversion '+logs[imin]
    print '... ls = '+repr(10**ls[imin])+' , lt = '+repr(10**lt[imin])
        
    
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
    from string import rjust
    
    
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
            P=load(datapath+rupt.split('.')[0]+'.'+rupt.split('.')[1]+'.sub'+rjust(str(k),4,'0')+'.stfpsd.npz')
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
                P=load(datapath+rupt.split('.')[0]+'.'+rupt.split('.')[1]+'.sub'+rjust(str(current[ks]),4,'0')+'.stfpsd.npz')
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
    from string import rjust
    
    dtopo=genfromtxt(dtopo_file)
    f=genfromtxt(fault_file)
    t=unique(dtopo[:,0])
    #Get z ranges
    z1=dtopo[:,3].min()
    z2=dtopo[:,3].max()
    zrange=max([-z1,z2])
    #Loop over time slices
    for kt in range(len(t)):
        print 't = '+str(t[kt])
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
        plt.savefig(out+rjust(str(kt),4,'0')+'vert.png')
        plt.close("all")
        

def plot_grd(grdfile,zlims,cmap,flip_lon=False):
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
    try:
        x=grd.variables['x'][:]
        y=grd.variables['y'][:]
        z=grd.variables['z'][:]
    except:
        x=grd.variables['lon'][:]
        y=grd.variables['lat'][:]
        z=grd.variables['z'][:]
    if flip_lon:
        x=x-360
    X,Y=meshgrid(x,y)
    plt.figure()
    plt.title(grdfile)
    plt.pcolormesh(X,Y,z,vmin=zlims[0],vmax=zlims[1],cmap=cmap)
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()




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