'''
Diego Melgar 01/2014
Green functions routines for source models
'''


        
def run_green(source,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax,smth):
    '''
    Compute GFs using Zhu & Rivera code for a given velocity model, source depth
    and station file. This function will make an external system call to fk.pl
    
    IN:
        source: 1-row numpy array containig informaiton aboutt he source, lat, lon, depth, etc...
        station_file: File name with the station coordinates
        dt: Desired sampling interval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon

    OUT:
        log: Sysytem standard output and standard error for log
    '''
    import subprocess
    from shlex import split
    
    depth='%.4f' % source[3]
    print("--> Computing GFs for source depth "+str(depth)+" km")
    #Get station distances to source
    d,az=src2sta(station_file,source)
    #Make distance string for system call
    diststr=''
    for k in range(len(d)):
        diststr=diststr+' %.3f' % d[k] #Truncate distance to 3 decimal palces (meters)
    if static==0: #Compute full waveform
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/'+str(smth)+'/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr)
        print("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/'+str(smth)+'/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr)
        print(command)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=p.communicate()
    else: #Compute only statics
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
        print("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
        print(command)
        p=subprocess.Popen(command,stdout=open('staticgf','w'),stderr=subprocess.PIPE)
        out,err=p.communicate()
    #Log output
    #print(out)
    #print(err)
    log=str(out)+str(err)
    return log
    
    

            
    
    
    
def okada_synthetics(strike,dip,rake,length,width,lon_source,lat_source,
                    depth_source,lon_obs,lat_obs,mu):
    '''
    Calculate neu synthetics for a subfault using Okada analytical solutions
    '''
    
    from okada_wrapper import dc3dwrapper
    from numpy import array,cos,sin,deg2rad,zeros
    from pyproj import Geod
    
    theta=strike-90
    theta=deg2rad(theta)
    
    #Rotaion matrices since okada_wrapper is only for east striking fault
    R=array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
    R2=array([[cos(-theta),-sin(-theta)],[sin(-theta),cos(-theta)]])
                       
    #position of point from lon/lat to x/y assuming subfault center is origin
    P=Geod(ellps='WGS84')
    az,baz,dist=P.inv(lon_source,lat_source,lon_obs,lat_obs)
    dist=dist/1000.
    x_obs=dist*sin(deg2rad(az))
    y_obs=dist*cos(deg2rad(az))
    
    #Calculate on rotated position
    xy=R.dot(array([x_obs, y_obs]))
    
    #Get Okada displacements
    lamb=mu
    alpha = (lamb + mu) / (lamb + 2 * mu)
    ss_in_m=1.0*cos(deg2rad(rake))
    ds_in_m=1.0*sin(deg2rad(rake))
    success, u, grad_u = dc3dwrapper(alpha, [xy[0], xy[1], 0.0],depth_source,dip,
                            [-length/2., length/2.], [-width/2., width/2],
                            [ss_in_m, ds_in_m, 0.0])
            
    #Rotate output
    urot=R2.dot(array([[u[0]], [u[1]]]))
    u[0]=urot[0]
    u[1]=urot[1]
    
    #output
    n=u[1]
    e=u[0]
    z=u[2]  
      
    return n,e,z
    
    

##########                   Utilities and stuff                      ##########          
        
def cartesian_azimuth(x,y,xs,ys):
    '''
    Compute source to station azimuths (from North) when sources given in cartesian coordinates
    
    IN:
        x: Vector of station x coordinates
        y: Vector of station y coordinates
        xs: Vector of source x coordinates
        ys: Vectpr of source y coordinates
        
    OUT:
        az: Source to station azimuth in degrees
    '''
    
    from numpy import arctan,nonzero,pi,intersect1d,rad2deg
    
    az=arctan((x-xs)/(y-ys))
    #Adjsut elements in 2nd quadrant
    ix=nonzero(x-xs<0)[0]
    iy=nonzero(y-ys>0)[0]
    i=intersect1d(ix,iy)
    az[i]=2*pi+az[i]
    #Adjust elements in 3rd and 4th quadrants
    i=nonzero(y-ys<0) 
    az[i]=pi+az[i]
    return rad2deg(az)
    
def src2sta(station_file,source,output_coordinates=False):
    '''
    Compute cartesian source to station distances and azimuths for all station/source pairs.
    
    IN:
        station_file: Path to station file
        source: numpy 1d array with source info read from file
        coord_type: =0 if coordinates are cartesian, =1 if they are lat/lon
    OUT:
        d - sorted distances vector in km
        az - azimuth from source to station in degrees
    '''
    
    from numpy import genfromtxt,zeros,array
    from obspy.geodetics.base import gps2dist_azimuth
    
    
    #Read station file
    #staname=genfromtxt(home+station_file,dtype="U",usecols=0)
    x=genfromtxt(station_file,dtype="f8",usecols=1)
    y=genfromtxt(station_file,dtype="f8",usecols=2)
    if x.shape==() or y.shape==(): #Single station file
        x=array([x])
        y=array([y])
    d=zeros(x.shape)
    az=zeros(x.shape)
    baz=zeros(x.shape)
    xs=source[1]
    ys=source[2]
    for k in range(len(x)):
        d[k],az[k],baz[k]=gps2dist_azimuth(ys,xs,y[k],x[k])
    d=d/1000
    
    if output_coordinates==True:
        return d,az,x,y
    else:
        return d,az
    
    

def origin_time(st,time_epi,tb,dt):
    '''
    Make start time of synthetics correspond with epicentral time
    
    Usage:
        st=origin_time(st,time_epi,tb)
    
    In:
        st: stream object to be altered
        time_epi: UTCDateTime object containing epicentral tiem
        tb: Number fo samples before first arrival in waveform
        
    Out:
        st: Time adjsuted waveform
    '''
    
    from datetime import timedelta

    p_wave_time = dt * tb + st[0].stats.sac.b #fk calculated time for p wave arrival time
    padding = dt * tb #padding time input
    t1 = time_epi + p_wave_time - padding #adjusted time for event, based on p wave arrival time and zero padding
    st[0].stats.starttime=t1
    #Set default sac headers to avoid invalid SAC write
    st[0].stats.sac['nzyear'] = t1.year
    st[0].stats.sac['nzjday'] = t1.julday
    st[0].stats.sac['nzhour'] = t1.hour
    st[0].stats.sac['nzmin'] = t1.minute
    st[0].stats.sac['nzsec'] = t1.second
    st[0].stats.sac['nzmsec'] = t1.microsecond/1000
    st[0].stats.sac.stla = 0.
    st[0].stats.sac.stlo = 0.
    return st

def rt2ne(r,t,azimuth):
    '''
    rotate time series of radial transverse to north and east. The azimuth is source 
    to station in degrees from north.
    '''
    from numpy import cos,sin,deg2rad
    az=deg2rad(azimuth)
    n=r*cos(az)-t*sin(az)
    e=r*sin(az)+t*cos(az)
    return n,e
    
def stdecimate(st,factor):
    '''
    Decimate stream by a constant factor, i.e. factor=4 will go from 4hz to 1Hz data
    '''
    from scipy.signal import filtfilt,butter
    
    #Anti-alias filter
    b, a = butter(10, 1./factor)
    y = filtfilt(b, a, st.data)
    stout=st.copy()
    stout.data=y
    #Decimate
    stout.decimate(factor,no_filter=True)
    return stout
    
def rtrim(st,T):
    '''
    Keep only the first T seconds of a waveform
    '''
    
    from datetime import timedelta
    
    stout=st.copy()
    start=st[0].stats.starttime
    T=timedelta(seconds=T)
    stout[0].trim(starttime=start,endtime=start+T)
    return stout
    
def triangle_stf(rise_time,dt):
    '''
    Make a triangle source time function of a given duration at a given sampling rate
    Area under triangle ahs to be one.
    '''
    
    from numpy import arange,r_
    
    rise_time=float(rise_time)
    t1=arange(0,rise_time/2,dt)
    t2=arange(rise_time/2,rise_time+dt,dt)
    
    m1=(4.*dt)/(rise_time**2)
    m2=-m1
    b2=(4.*dt)/rise_time
    
    y1=m1*t1
    y2=m2*t2+b2
    
    t=r_[t1,t2]
    stf=r_[y1,y2]
    return t,stf
    
def dreger_stf(rise_time,zeta,dt):
    '''
    '''
    from numpy import arange
    
    rise_time=float(rise_time)
    t=arange(0,rise_time/2,dt)
    
    

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured

    
    
    

