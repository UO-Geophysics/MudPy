'''
Diego Melgar 01/2014
Green functions routines for source models
'''


        
def run_green(source,station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax):
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
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/1/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr)
        print("fk.pl -M"+model_name+"/"+depth+"/f -N"+str(NFFT)+"/"+str(dt)+'/1/'+repr(dk)+' -P'+repr(pmin)+'/'+repr(pmax)+'/'+repr(kmax)+diststr)
        #print(command)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=p.communicate()
    else: #Compute only statics
        command=split("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
        print("fk.pl -M"+model_name+"/"+depth+"/f -N1 "+diststr)
        #print(command)
        p=subprocess.Popen(command,stdout=open('staticgf','w'),stderr=subprocess.PIPE)
        out,err=p.communicate()
    #Log output
    #print(out)
    #print(err)
    log=str(out)+str(err)
    return log
    
    

            
    
    
    
def run_syn(home,project_name,source,station_file,green_path,model_name,integrate,static,tsunami,
        subfault,time_epi,beta,impulse=False,okada=False,okada_mu=45e9,insar=False):
    '''
    Use green functions and compute synthetics at stations for a single source
    and multiple stations. This code makes an external system call to syn.c first it
    will make the external call for the strike-slip component then a second externall
    call will be made for the dip-slip component. The unit amount of moment is 1e15
    which corresponds to Mw=3.9333...
    
    IN:
        source: 1-row numpy array containig informaiton aboutt he source, lat, lon, depth, etc...
        station_file: File name with the station coordinates
        green_path: Directopry where GFs are stored
        model_file: File containing the Earth velocity structure
        integrate: =0 if youw ant velocity waveforms, =1 if you want displacements
        static: =0 if computing full waveforms, =1 if computing only the static field
        subfault: String indicating the subfault being worked on
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        log: Sysytem standard output and standard error for log
    '''

    import os
    import subprocess
    from mudpy.forward import get_mu
    from numpy import array,genfromtxt,loadtxt,savetxt,log10,argmin
    from obspy import read
    from shlex import split
    
    #Constant parameters
    rakeDS=90+beta #90 is thrust, -90 is normal
    rakeSS=0+beta #0 is left lateral, 180 is right lateral
    tb=50 #Number of samples before first arrival
    #Load structure
    model_file=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_file,ndmin=2)
    #Parse the soruce information
    num=str(int(source[0])).rjust(4,'0')
    xs=source[1]
    ys=source[2]
    zs=source[3]
    strike=source[4]
    dip=source[5]
    rise=source[6]
    if impulse==True:  #Impulse GFs or apply rise time
        duration=0
    else:
        duration=source[7]
    ss_length=source[8]
    ds_length=source[9]
    ss_length_in_km=ss_length/1000.
    ds_length_in_km=ds_length/1000.
    strdepth='%.4f' % zs
    if static==0 and tsunami==0:  #Where to save dynamic waveforms
        green_path=green_path+'dynamic/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    if static==0 and tsunami==1:  #Where to save dynamic waveforms
        green_path=green_path+'tsunami/'+model_name+"_"+strdepth+".sub"+subfault+"/"
    print("--> Computing synthetics at stations for the source at ("+str(xs)+" , "+str(ys)+")")
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    if staname.shape==(): #Single staiton file
        staname=array([staname])
    #Compute distances and azimuths
    d,az,lon_sta,lat_sta=src2sta(station_file,source,output_coordinates=True)
    #Get moment corresponding to 1 meter of slip on subfault
    mu=get_mu(structure,zs)
    Mo=mu*ss_length*ds_length*1
    Mw=(2./3)*(log10(Mo)-9.1)
    #Move to output folder
    log='' #Initalize log
    os.chdir(green_path)
    for k in range(len(d)):
        if static==0: #Compute full waveforms
            diststr='%.3f' % d[k] #Need current distance in string form for external call
            #Form the strings to be used for the system calls according to suer desired options
            if integrate==1: #Make displ.
                #First Stike-Slip GFs
                commandSS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.disp.x -G"+green_path+diststr+".grn.0"
                print(commandSS) #Output to screen so I know we're underway
                log=log+commandSS+'\n' #Append to log
                commandSS=split(commandSS) #Split string into lexical components for system call
                #Now dip slip
                commandDS="syn -I -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.disp.x -G"+green_path+diststr+".grn.0"
                print(commandDS)
                log=log+commandDS+'\n'
                commandDS=split(commandDS)
            else: #Make vel.
                #First Stike-Slip GFs
                commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".SS.vel.x -G"+green_path+diststr+".grn.0"
                print(commandSS)
                log=log+commandSS+'\n'
                commandSS=split(commandSS)
                #Now dip slip
                commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+" -D"+str(duration)+ \
                    "/"+str(rise)+" -A"+str(az[k])+" -O"+staname[k]+".subfault"+num+".DS.vel.x -G"+green_path+diststr+".grn.0"
                print(commandDS)
                log=log+commandDS+'\n'
                commandDS=split(commandDS)
            #Run the strike- and dip-slip commands (make system calls)
            p=subprocess.Popen(commandSS,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate() 
            p=subprocess.Popen(commandDS,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            out,err=p.communicate()
            #Result is in RTZ system (+Z is down) rotate to NEZ with +Z up and scale to m or m/s
            if integrate==1: #'tis displacememnt
                #Strike slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.SS.disp.r')
                    t=read(staname[k]+".subfault"+num+'.SS.disp.t')
                    z=read(staname[k]+".subfault"+num+'.SS.disp.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.SS.disp.ri')
                    t=read(staname[k]+".subfault"+num+'.SS.disp.ti')
                    z=read(staname[k]+".subfault"+num+'.SS.disp.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                #Scale to m and overwrite with rotated waveforms
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.SS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.disp.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.SS.disp.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.SS.disp.r')
                silentremove(staname[k]+".subfault"+num+'.SS.disp.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.ri')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.ti')
                    silentremove(staname[k]+".subfault"+num+'.SS.disp.zi')
                #Dip Slip
                if duration>0:
                    r=read(staname[k]+".subfault"+num+'.DS.disp.r')
                    t=read(staname[k]+".subfault"+num+'.DS.disp.t')
                    z=read(staname[k]+".subfault"+num+'.DS.disp.z')
                else:
                    r=read(staname[k]+".subfault"+num+'.DS.disp.ri')
                    t=read(staname[k]+".subfault"+num+'.DS.disp.ti')
                    z=read(staname[k]+".subfault"+num+'.DS.disp.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.DS.disp.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.disp.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.DS.disp.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.DS.disp.r')
                silentremove(staname[k]+".subfault"+num+'.DS.disp.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.ri')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.ti')
                    silentremove(staname[k]+".subfault"+num+'.DS.disp.zi')
            else: #Waveforms are velocity, as before, rotate from RT-Z to NE+Z and scale to m/s
                #Strike slip
                if duration>0: #Is there a source time fucntion? Yes!
                    r=read(staname[k]+".subfault"+num+'.SS.vel.r')
                    t=read(staname[k]+".subfault"+num+'.SS.vel.t')
                    z=read(staname[k]+".subfault"+num+'.SS.vel.z')
                else: #No! This is the impulse response!
                    r=read(staname[k]+".subfault"+num+'.SS.vel.ri')
                    t=read(staname[k]+".subfault"+num+'.SS.vel.ti')
                    z=read(staname[k]+".subfault"+num+'.SS.vel.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.SS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.SS.vel.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.SS.vel.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.SS.vel.r')
                silentremove(staname[k]+".subfault"+num+'.SS.vel.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.ri')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.ti')
                    silentremove(staname[k]+".subfault"+num+'.SS.vel.zi')
                #Dip Slip
                if duration>0:
                    r=read(staname[k]+".subfault"+num+'.DS.vel.r')
                    t=read(staname[k]+".subfault"+num+'.DS.vel.t')
                    z=read(staname[k]+".subfault"+num+'.DS.vel.z')
                else:
                    r=read(staname[k]+".subfault"+num+'.DS.vel.ri')
                    t=read(staname[k]+".subfault"+num+'.DS.vel.ti')
                    z=read(staname[k]+".subfault"+num+'.DS.vel.zi')
                ntemp,etemp=rt2ne(r[0].data,t[0].data,az[k])
                n=r.copy()
                n[0].data=ntemp/100
                e=t.copy()
                e[0].data=etemp/100
                z[0].data=z[0].data/100
                n=origin_time(n,time_epi,tb)
                e=origin_time(e,time_epi,tb)
                z=origin_time(z,time_epi,tb)
                n.write(staname[k]+".subfault"+num+'.DS.vel.n',format='SAC')
                e.write(staname[k]+".subfault"+num+'.DS.vel.e',format='SAC')
                z.write(staname[k]+".subfault"+num+'.DS.vel.z',format='SAC')
                silentremove(staname[k]+".subfault"+num+'.DS.vel.r')
                silentremove(staname[k]+".subfault"+num+'.DS.vel.t')
                if impulse==True:
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.ri')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.ti')
                    silentremove(staname[k]+".subfault"+num+'.DS.vel.zi')
        else: #Compute static synthetics
            os.chdir(green_path+'static/') #Move to appropriate dir
            if okada==False:
                diststr='%.1f' % d[k] #Need current distance in string form for external call
                
                
                if insar==True:
                    green_file=model_name+".static."+strdepth+".sub"+subfault+'.insar' #Output dir
                else: #GPS
                    green_file=model_name+".static."+strdepth+".sub"+subfault+'.gps' #Output dir
                
                
                print(green_file)
                log=log+green_file+'\n' #Append to log
                statics=loadtxt(green_file) #Load GFs
                #Print static GFs into a pipe and pass into synthetics command
                station_index=argmin(abs(statics[:,0]-d[k])) #Look up by distance
                try:
                    temp_pipe=statics[station_index,:]
                except:
                    temp_pipe=statics
                inpipe=''
                for j in range(len(temp_pipe)):
                    inpipe=inpipe+' %.6e' % temp_pipe[j]
                #Form command for external call
                commandDS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeDS)+\
                        " -A"+str(az[k])+" -P"
                commandSS="syn -M"+str(Mw)+"/"+str(strike)+"/"+str(dip)+"/"+str(rakeSS)+\
                        " -A"+str(az[k])+" -P"
                print(staname[k])
                print(commandSS)
                print(commandDS)
                log=log+staname[k]+'\n'+commandSS+'\n'+commandDS+'\n' #Append to log
                commandSS=split(commandSS) #Lexical split
                commandDS=split(commandDS)
                #Make system calls, one for DS, one for SS, and save log
                ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
                p=subprocess.Popen(commandSS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.SS.static.rtz','w'),stderr=subprocess.PIPE)     
                out,err=p.communicate()  
                log=log+str(out)+str(err)
                ps=subprocess.Popen(['printf',inpipe],stdout=subprocess.PIPE,stderr=subprocess.PIPE)  #This is the statics pipe, pint stdout to syn's stdin
                p=subprocess.Popen(commandDS,stdin=ps.stdout,stdout=open(staname[k]+'.subfault'+num+'.DS.static.rtz','w'),stderr=subprocess.PIPE)     
                out,err=p.communicate() 
                log=log+str(out)+str(err)       
                #Rotate radial/transverse to East/North, correct vertical and scale to m
                statics=loadtxt(staname[k]+'.subfault'+num+'.SS.static.rtz')
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
                statics=loadtxt(staname[k]+'.subfault'+num+'.DS.static.rtz')
                u=statics[2]/100
                r=statics[3]/100
                t=statics[4]/100
                ntemp,etemp=rt2ne(array([r,r]),array([t,t]),az[k])
                n=ntemp[0]
                e=etemp[0]
                savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
            else:
                #SS
                n,e,u=okada_synthetics(strike,dip,rakeSS,ss_length_in_km,ds_length_in_km,xs,ys,
                    zs,lon_sta[k],lat_sta[k],okada_mu)
                savetxt(staname[k]+'.subfault'+num+'.SS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
                #DS
                n,e,u=okada_synthetics(strike,dip,rakeDS,ss_length_in_km,ds_length_in_km,xs,ys,
                    zs,lon_sta[k],lat_sta[k],okada_mu)
                savetxt(staname[k]+'.subfault'+num+'.DS.static.neu',(n,e,u,beta),header='north(m),east(m),up(m),beta(degs)')
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
    
    

def origin_time(st,time_epi,tb):
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
    
    t1=st[0].stats.starttime  #Waveform starttime
    td=timedelta(seconds=st[0].stats.delta*tb)  #Shift due to pre-first arrival samples
    #Shift forward
    t1=t1+td
    #Shift to oring time
    t1=time_epi+timedelta(minutes=t1.minute,seconds=t1.second,microseconds=t1.microsecond)-td
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

    
    
    

