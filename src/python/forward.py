'''
D. Melgar 02/2014

Forward modeling routines
'''

def waveforms(home,project_name,rupture_name,station_file,model_name,integrate):
    '''
    '''
    from numpy import loadtxt,genfromtxt,deg2rad,sin,cos,allclose
    from obspy import read,Stream
    from string import rjust
    import datetime
    
    #constants
    unitM=1e15 #N-m
    
    #Output where?
    outpath=home+project_name+'/output/forward_models/'
    logpath=home+project_name+'/logs/'
    log=''
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #load source
    source=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    #Load stations
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    #Load structure
    model_path=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_path)
    #What am I processing v or d?
    if integrate==1:
        vord='disp'
    else:
        vord='vel'
        
    #Loop over stations
    for ksta in range(len(staname)):
        print 'Working on station '+staname[ksta]+' ('+str(ksta)+'/'+str(len(staname))+')'
        #Initalize output
        n=Stream()
        e=Stream()
        z=Stream()
        sta=staname[ksta]
        #Loop over sources (Add delays)
        for k in range(source.shape[0]):
            #Get subfault parameters
            nfault='subfault'+rjust(str(int(source[k,0])),4,'0')
            nsub='sub'+rjust(str(int(source[k,0])),4,'0')
            zs=source[k,3]
            rake=source[k,6]
            slip=source[k,9]
            sslength=source[k,10]
            dslength=source[k,11]
            rtime=source[k,12]
            #Where's the data
            strdepth='%.4f' % zs 
            syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
            #Find rigidity at this source point
            mu=get_mu(structure,zs)
            #Compute equivalent Moment at this source point
            Mo=mu*slip*sslength*dslength
            #Get synthetics
            ess=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.e')
            nss=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.n')
            zss=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.z')
            eds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.e')
            nds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.n')
            zds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.z')
            #Time shift them according to subfault rupture time
            dt=ess[0].stats.delta
            ess=tshift(ess,rtime)
            ess[0].stats.starttime=round_time(ess[0].stats.starttime,dt)
            nss=tshift(nss,rtime)
            nss[0].stats.starttime=round_time(nss[0].stats.starttime,dt)
            zss=tshift(zss,rtime)
            zss[0].stats.starttime=round_time(zss[0].stats.starttime,dt)
            eds=tshift(eds,rtime)
            eds[0].stats.starttime=round_time(eds[0].stats.starttime,dt)
            nds=tshift(nds,rtime)
            nds[0].stats.starttime=round_time(nds[0].stats.starttime,dt)
            zds=tshift(zds,rtime)
            zds[0].stats.starttime=round_time(zds[0].stats.starttime,dt)
            #get rake contribution and moment multiplier
            dsmult=sin(deg2rad(rake))
            ssmult=cos(deg2rad(rake))
            M=Mo/unitM
            if allclose(slip,0)==False:  #Only add things that matter
                log=log+nfault+', SS='+str(ssmult)+', DS='+str(dsmult)+', Mscale='+str(M)+'\n'
                #A'ight, add 'em up
                etotal=add_traces(ess,eds,ssmult,dsmult,M)
                ntotal=add_traces(nss,nds,ssmult,dsmult,M)
                ztotal=add_traces(zss,zds,ssmult,dsmult,M)
                #Add to previous subfault's results
                e=add_traces(e,etotal,1,1,1)
                n=add_traces(n,ntotal,1,1,1)
                z=add_traces(z,ztotal,1,1,1)
            else:
                log=log+"No slip on subfault "+nfault+', ignoring it...\n'
                
        #Save results
        e[0].data=e[0].data.filled()  #This is a workaround of a bug in obspy
        n[0].data=n[0].data.filled()
        z[0].data=z[0].data.filled()
        e.write(outpath+sta+'.'+vord+'.e',format='SAC')
        n.write(outpath+sta+'.'+vord+'.n',format='SAC')
        z.write(outpath+sta+'.'+vord+'.z',format='SAC')
    f=open(logpath+'waveforms.'+now+'.log','a')
    f.write(log)
    f.close()
        

def coseismics(home,project_name,rupture_name,station_file,model_name):
    '''
    '''
    from numpy import loadtxt,genfromtxt,deg2rad,sin,cos,array,savetxt
    from string import rjust
    
    #constants
    unitM=1e15 #N-m
    
    #Output where?
    outpath=home+project_name+'/output/forward_models/'
    #load source
    source=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    #Load stations
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    #Load structure
    model_path=home+project_name+'/structure/'+model_name
    structure=loadtxt(model_path)
    #Loop over stations
    for ksta in range(len(staname)):
        #Initalize output
        n=array([0])
        e=array([0])
        z=array([0])
        sta=staname[ksta]
        print 'Working on station'+sta
        #Loop over sources
        for k in range(source.shape[0]):
            #Get subfault parameters
            nfault='subfault'+rjust(str(int(source[k,0])),4,'0')
            zs=source[k,3]
            rake=source[k,6]
            slip=source[k,9]
            sslength=source[k,10]
            dslength=source[k,11]
            #Where's the data
            syn_path=home+project_name+'/GFs/static/'
            #Find rigidity at this source point
            mu=get_mu(structure,zs)
            #Compute equivalent Moment at this source point
            Mo=mu*slip*sslength*dslength
            #Get synthetics
            coseis_ss=loadtxt(syn_path+sta+'.'+nfault+'.SS.static.enu')
            ess=coseis_ss[0]
            nss=coseis_ss[1]
            zss=coseis_ss[2]
            coseis_ds=loadtxt(syn_path+sta+'.'+nfault+'.DS.static.enu')
            eds=coseis_ds[0]
            nds=coseis_ds[1]
            zds=coseis_ds[2]
            #get rake contribution and moment multiplier
            dsmult=sin(deg2rad(rake))
            ssmult=cos(deg2rad(rake))
            M=Mo/unitM
            #A'ight, add 'em up
            etotal=M*(dsmult*eds+ssmult*ess)
            ntotal=M*(dsmult*nds+ssmult*nss)
            ztotal=M*(dsmult*zds+ssmult*zss)
            #Add to previous subfault's results
            e=e+etotal
            n=n+ntotal
            z=z+ztotal
        #Save results
        savetxt(outpath+sta+'.static.enu',(e,n,z))
            
            
            
###########                Tools and trinkets                      #############
    
def get_mu(structure,zs):
    '''
    Look in structure model and compute rigidity
    '''
    from numpy import nonzero
    
    Z=structure[:,0].cumsum()
    #Which layer do I want?
    i=nonzero(zs>Z)[0]
    if i.size==0: #It's in top layer
        imu=0
    else:
        imu=max(i)+1
    if imu>=structure.shape[0]:
        imu=imu-1#It's int he half-space
    mu=((1000*structure[imu,1])**2)*structure[imu,3]
    #print "Rigidity at z="+str(zs)+' is, mu = '+str(mu/1e9)+'GPa'
    return mu
    
def add_traces(ss,ds,ssmult,dsmult,M):
    '''
    Add two stream objects with dip slip and striek slip contributions
    
    For simple addition use ss=ds=M=1
    
    NOTES: Right now I'm truncating tiems tot eh nearest second, this si incorrect, they shoudl be truncated tot he enarest dt interval.
    Also pad value should be zero for VELOCITY and the last sampelf or DISPALCEMMENT, this needs adjsutment
    '''
    from numpy import zeros
    
    #If one stream object is empty set it to zeros, if bot are empty then freak out
    if ss.count()==0 and ds.count()==0:
        print 'FATAL ERROR: can\'t add 2 empty stream objects doofus'
        return None
    if ss.count()==0:
        ss=ds.copy()
        ss[0].data=zeros(ds[0].data.shape)
    elif ds.count()==0:
        ds=ss.copy()
        ds[0].data=zeros(ss[0].data.shape) 
        
    #Round times to dt sampling interval
    ss[0].stats.starttime=round_time(ss[0].stats.starttime,ss[0].stats.delta)
    ds[0].stats.starttime=round_time(ds[0].stats.starttime,ds[0].stats.delta)
    #Find earliest starttime
    if ss[0].stats.starttime<=ds[0].stats.starttime:
        t1=ss[0].stats.starttime
    else:
        t1=ds[0].stats.starttime
    #Find altest end time
    if ss[0].stats.endtime>=ds[0].stats.endtime:
        t2=ss[0].stats.endtime
    else:
        t2=ds[0].stats.endtime
    #Now extend both arrays and fill start with zeros then end with last sample
    ss[0].trim(t1,ss[0].stats.endtime,pad=True,fill_value=0)
    ds[0].trim(t1,ds[0].stats.endtime,pad=True,fill_value=0)
    fillend=ss[0].data[-60:-1].mean()
    ss[0].trim(t1,t2,pad=True,fill_value=fillend)
    fillend=ds[0].data[-20:-1].mean()
    ds[0].trim(t1,t2,pad=True,fill_value=fillend)
    #Apply rake scaling value
    ss[0].data=ss[0].data*ssmult
    ds[0].data=ds[0].data*dsmult
    #Creat output stream
    st=ss.copy()
    #Add and apply scale
    st[0].data=M*(ss[0].data+ds[0].data)
    #And done
    return st


def tshift(st,tshift):
    '''
    shift a stream object by tshift seconds, positive moves forward in time
    '''
    from datetime import timedelta
    td=timedelta(seconds=tshift)
    st[0].stats.starttime=st[0].stats.starttime+td
    return st
        
        
def round_time(t1,delta):
    '''
    '''
    from datetime import timedelta
    #Move start and end times to start exactly on the dt intervals
    dtmicro=delta*1e6
    intervals=t1.microsecond/dtmicro #How many intervals in microseconds
    adjustment=(round(intervals)-intervals)*dtmicro
    td=timedelta(microseconds=adjustment)
    t1=t1+td
    return t1
