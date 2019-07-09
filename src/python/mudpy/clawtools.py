def make_dtopo():
    '''
    Make geoclaw dtopo file
    '''
    from numpy import genfromtxt,zeros
    
    #Run params
    f='/Users/dmelgarm/Research/Slip_Inv/tohoku_tsunami/'
    stafile='tohoku.sta'
    dlon=0.033333
    dlat=0.033333
    dt=5
    stat_or_dyn='s'
    
    #Get station list
    
    sta=genfromtxt(f+'data/station_info/'+stafile,usecols=0,dtype='S4')
    s=genfromtxt(f+'data/station_info/'+stafile,usecols=[1,2])
    lon=s[:,0]
    lat=s[:,1]
    if stat_or_dyn.lower()=='s':
        n=zeros(len(sta))
        e=n.copy()
        u=n.copy()
        for ksta in range(len(sta)):
            print(ksta)
            neu=genfromtxt(f+'output/forward_models/'+str(sta[ksta])+'.static.neu')
            n[ksta]=neu[0]
            e[ksta]=neu[1]
            u[ksta]=neu[2]
            print(neu[2])

def make_grid(lon1,lon2,lat1,lat2,dlon,dlat,prefix,outfile):
    '''
    Make a grid of points where you will request coseismic deformation
    
    lon1,lon2,lat1,lat2 are grid limits
    dlon,dlat are lon and lat spacing
    prefix is 2 character "station name" to give seafloor points
    outfile is absolute path to file name
    '''
    
    from numpy import arange
    
    lon=arange(lon1,lon2+dlon,dlon)
    lat=arange(lat1,lat2+dlat,dlat)
    lat=lat[::-1]
    #Turn into pseudo station file
    sta=arange(0,len(lon)*len(lat),1)
    ksta=0
    fid=open(outfile,'w')
    for klat in range(len(lat)):
        for klon in range(len(lon)):
            outsta=prefix+str(sta[ksta]).rjust(4,'0')
            outlon='%.4f' % lon[klon]
            outlat='%.4f' % lat[klat]
            ksta+=1
            out=outsta+'\t'+outlon+'\t'+outlat+'\n'
            fid.write(out)
    fid.close()    
    
#Plot a transect of eta.blend files
def plot_transect(path,N,tsunlims,dlims,dt,xyannot):
    '''
    Plot snapshots along a transect
    '''
    from obspy.core.util.geodetics import gps2DistAzimuth
    from numpy import genfromtxt,zeros
    from matplotlib import pyplot as plt
    import matplotlib
    
    font = {'family' : 'normal','size'   : 14}

    matplotlib.rc('font', **font)
    trench_index=395 #Where is trench
    cut_index=160 #Where is land
    trench_index=trench_index-cut_index
    fig, axarr = plt.subplots(N, 1) 
    for k in range(N):
        time=k*dt
        trackname='AA.f'+str(k+1).rjust(4,'0')+'.track'
        t=genfromtxt(path+trackname)
        t=t[cut_index:,:]
        lon1=t[trench_index,0]
        lat1=t[trench_index,1]
        d=zeros(len(t))
        for n in range(len(t)):
            d[n],Az,BAz=gps2DistAzimuth(lat1,lon1,t[n,1],t[n,0])
        d=d/1000
        d[0:trench_index]=-d[0:trench_index]
        tstring=str(time)+'s'
        #Plot
        ax=axarr[k]
        ax.plot(d,t[:,2],'r',lw=2)
        ax.set_ylim(tsunlims)
        ax.set_xlim(dlims)
        ax.grid(which='both')
        ax.yaxis.set_ticklabels(['',0,'',10,'',20,''])
        if k!=N-1:
            ax.xaxis.set_ticklabels([])
        ax.annotate(tstring,xy=xyannot,fontsize=14)
        if k==N-1:
            ax.xaxis.set_label('Distance from trench (km)')
            print(t[0,0])
            print(t[-1,0])
            print(t[0,1])
            print(t[-1,1])
    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.8, top=0.9, wspace=0, hspace=0)
        
       
def gauge2sac(gauge_file,dictionary,xyfile,outdir,time_epi,dt):
    '''
    Convert output from fort.gauge file intop individual sac files
    ''' 
    from numpy import genfromtxt,unique,where,arange,interp
    from obspy import Stream,Trace
    from obspy.core.util.attribdict import AttribDict

    
    #Read gauge file
    gauges=genfromtxt(gauge_file)
    #Read names
    plume_name=genfromtxt(dictionary,usecols=0,dtype='S')
    claw_name=genfromtxt(dictionary,usecols=1)
    lat=genfromtxt(xyfile,usecols=2)
    lon=genfromtxt(xyfile,usecols=3)
    #Find unique stations
    gauge_list=unique(gauges[:,0])
    for k in range(len(gauge_list)):
        print(k)
        st=Stream(Trace())
        i=where(gauges[:,0]==gauge_list[k])[0]
        data=gauges[i,6]
        time=gauges[i,2]
        ti=arange(0,time.max(),dt)
        tsunami=interp(ti,time,data)
        st[0].data=tsunami
        st[0].stats.starttime=time_epi
        st[0].stats.delta=dt
        iname=where(claw_name==gauge_list[k])[0][0]
        st[0].stats.station=plume_name[iname]
        sac=AttribDict()
        sac.stla=lat[iname]
        sac.stlo=lon[iname]
        sac.evla=46.607
        sac.evlo=153.230
        #sac.iztype='IO'
        st[0].stats['sac']=sac
        st.write(outdir+'/'+plume_name[iname]+'.tsun.sac',format='SAC')
        

        
        
        
    
    
    