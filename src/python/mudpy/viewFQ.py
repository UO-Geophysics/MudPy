'''
Module for plotting FakeQuakes related studd
'''

def plot_LW_scaling(home,project_name,run_name):
    '''
    Plot a comparisson between used length/width and Blasser scaling laws
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10
    from matplotlib import pyplot as plt
    
    logs=glob(home+project_name+'/output/ruptures/'+run_name+'*.log')
    L=zeros(len(logs))
    W=zeros(len(logs))
    Mw=zeros(len(logs))
    for k in range(len(logs)):
        f=open(logs[k],'r')
        while True:
            line=f.readline()
            if 'Lmax' in line:
                L[k]=float(line.split(':')[-1].split()[0])
            if 'Wmax' in line:
                W[k]=float(line.split(':')[-1].split()[0])
            if 'Actual magnitude' in line:
                Mw[k]=float(line.split()[-1])
                break
    
    #Make plot
    plt.figure(figsize=(9,4))
    
    plt.subplot(121)
    Mw_synth=arange(7.8,9.3,0.1)
    plt.plot(Mw_synth,-2.37+0.57*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(L),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (L) [km]')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(8.0,3.17))
    
    plt.subplot(122)
    plt.plot(Mw_synth,-1.86+0.46*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(W),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (W) [km]')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(8.0,2.8))
    plt.show()
    
    plt.subplots_adjust(bottom=0.15)
    

def one_event_pgd_scaling(summary_file,event_log):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,log10
    from matplotlib import pyplot as plt
    from string import replace
    
    # Read summary file
    lonlat=genfromtxt(summary_file,usecols=[1,2])
    pgd=genfromtxt(summary_file,usecols=[6])*100
    
    # Get hypocenter
    f=open(event_log,'r')
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

    #compute station to hypo distances
    d=zeros(len(lonlat))
    for k in range(len(lonlat)):
        d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
        d[k]=d[k]/1000
        
    fig = plt.figure()
    ax = plt.gca()
    ax.scatter(d,pgd)
    ax.set_yscale('log')
    ax.set_xscale('log')
    # Plot reference line
    dref=logspace(0,3,50)
    pgdref=10**(-4.434+1.047*Mw-0.138*Mw*log10(dref))
    plt.plot(dref,pgdref)
    plt.show()
    
    
def record_section(home,project_name,GF_list,rupture,factor=10):
    '''
    Plot record section
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,ones,zeros
    from matplotlib import pyplot as plt
    from string import replace
    from obspy import read
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    xmajorLocator = MultipleLocator(50)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(200)
    yminorLocator = MultipleLocator(50)
    
    # Read summary file
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    event_log=home+project_name+'/output/ruptures/'+rupture+'.log'
    
    #Load velocity model for ray tracing
    velmod = velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
    
    # Get hypocenter
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'Hypocenter (lon,lat,z[km])' in line:
            s=replace(line.split(':')[-1],'(','')
            s=replace(s,')','')
            hypo=array(s.split(',')).astype('float')
            loop_go=False

    #compute station to hypo distances
    d=zeros(len(lonlat))
    for k in range(len(lonlat)):
        d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
        d[k]=d[k]/1000
        if lonlat[k,1]<hypo[1]: #station is south
            d[k]=-d[k]
 
    ptime=1e6*ones(len(sta))
    stime=1e6*ones(len(sta))
    for k in range(len(sta)):
        offset=d[k]
        e=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYE.sac')
        n=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYN.sac')
        z=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYZ.sac')            
        
        #Normalize and scale and apply distance offset
        n[0].data=factor*(n[0].data/max(abs(n[0].data)))+offset
        e[0].data=factor*(e[0].data/max(abs(e[0].data)))+offset
        z[0].data=factor*(z[0].data/max(abs(z[0].data)))+offset
        
        if k==0:
            N=n.copy()
            E=e.copy()
            Z=z.copy()
        else:
            N+=n
            E+=e
            Z+=z
        
        #Get ray arrivals
        deg=locations2degrees(hypo[1],hypo[0],lonlat[k,1],lonlat[k,0])
        
        # Ray trace
        try: #This is to deal with some weird bug in obspy.TauP
            arrivals = velmod.get_travel_times(source_depth_in_km=hypo[2],distance_in_degree=deg,phase_list=['P','Pn','S','Sn','p','s'])
        except:
            print 'No phases for '+rupture
            plt.close()
            return
        
        #Determine P and S arrivals
        for kphase in range(len(arrivals)):
            if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
                if arrivals[kphase].time<ptime[k]:
                    ptime[k]=arrivals[kphase].time
            if 'S' == arrivals[kphase].name or 's' == arrivals[kphase].name or 'Sn' == arrivals[kphase].name:
                if arrivals[kphase].time<stime[k]:
                    stime[k]=arrivals[kphase].time
    
    #plot on figure
    plt.figure(figsize=(18,7))
    
    plt.subplot(131)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(N[k].times(),N[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([n[0].times()[0],n[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    plt.ylabel('Hypocentral distance (km)')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.title('North')
    
    plt.subplot(132)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(E[k].times(),E[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([e[0].times()[0],e[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.tick_params(axis='y',labelleft='off')
    plt.title('East')
    
    plt.subplot(133)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(Z[k].times(),Z[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([z[0].times()[0],z[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.tick_params(axis='y',labelleft='off')
    plt.title('Vertical')
    
    plt.subplots_adjust(left=0.05,right=0.98,wspace=0.05)
    plt.show()
    
    
    
def statics_map(summary_file,scale,xl,yl,vertscale=400):
    '''
    Plot statics
    '''
    
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,zeros,r_,where
    
    sta=genfromtxt(summary_file)
    
    
    plt.figure(figsize=(5,12))
    
    #Plot coasts (also hard coded to cascadia for now)
    fcoast=open('/Users/dmelgar/Cascadia/coast/coast.txt','r')
    line=fcoast.readline()
    parsing_coast=True
    while parsing_coast==True:
        if '>' in line:
            complete_polygon=False
            count=0
            while complete_polygon==False:
                line=fcoast.readline()
                if '>' in line or line=='':
                    complete_polygon=True
                    plt.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.3,c='k')
                else:
                    new_coordinates=zeros((1,2))
                    new_coordinates[0,0]=float(line.split()[0])
                    new_coordinates[0,1]=float(line.split()[1])
                    if count==0:
                        coast_coordinates=new_coordinates.copy()
                        count+=1
                    else:
                        coast_coordinates=r_[coast_coordinates,new_coordinates]
                if line=='':
                    parsing_coast=False
    #Done plotting coast
    #Plot cnational boundaries(also hard coded to cascadia for now)
    fcoast=open('/Users/dmelgar/Cascadia/coast/boundaries.txt','r')
    line=fcoast.readline()
    parsing_coast=True
    while parsing_coast==True:
        if '>' in line:
            complete_polygon=False
            count=0
            while complete_polygon==False:
                line=fcoast.readline()
                if '>' in line or line=='':
                    complete_polygon=True
                    plt.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.5,c='k')
                else:
                    new_coordinates=zeros((1,2))
                    new_coordinates[0,0]=float(line.split()[0])
                    new_coordinates[0,1]=float(line.split()[1])
                    if count==0:
                        coast_coordinates=new_coordinates.copy()
                        count+=1
                    else:
                        coast_coordinates=r_[coast_coordinates,new_coordinates]
                if line=='':
                    parsing_coast=False
    #Done plotting coast
    
    iup=where(sta[:,5]>0)
    idown=where(sta[:,5]<0)
    plt.scatter(sta[iup,1],sta[iup,2],c='r',s=sta[iup,5]*vertscale,lw=0.1)
    plt.scatter(sta[idown,1],sta[idown,2],c='#4169E1',s=-sta[idown,5]*vertscale,lw=0.1)
    plt.quiver(sta[:,1],sta[:,2],sta[:,4],sta[:,3],scale=scale,width=0.007)
    
    plt.xlim(xl)
    plt.ylim(yl)
    
    #Reference
    ref=1
    plt.quiver(-123.8,50.7,ref,0,scale=scale,width=0.007)
    plt.annotate('1m',xy=(-121.5,50.6))
    
    ref=1
    plt.scatter(-122.6,50.1,c='r',s=ref*vertscale,lw=0.1)
    plt.quiver(-122.6,49.98,0,1,10,width=0.007)
    plt.scatter(-123.6,50.1,c='#4169E1',s=ref*vertscale,lw=0.1)
    plt.quiver(-123.6,50.25,0,-1,10,width=0.007)
    plt.annotate('1m',xy=(-121.5,50.0))
    
    plt.show()
    


    

    