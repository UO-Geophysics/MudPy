'''
D. Melgar 02/2014

Forward modeling routines
'''

def waveforms(home,project_name,rupture_name,station_file,model_name,integrate):
    '''
    '''
    from numpy import loadtxt,genfromtxt,deg2rad,sin,cos,round,allclose
    from obspy import read,Stream
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
    #What am I processing v or d
    if integrate==1:
        vord='disp'
    else:
        vord='vel'
        
    #Loop over stations
    for ksta in range(len(staname)):
        #Initalize output
        n=Stream()
        e=Stream()
        z=Stream()
        sta=staname[ksta]
        #Loop over sources (Add delays)
        for k in range(source.shape[0]):
            #Get subfault parameters
            nfault='subfault'+rjust(str(int(source[k,0])),4,'0')
            zs=source[k,3]
            rake=source[k,6]
            slip=source[k,9]
            sslength=source[k,10]
            dslength=source[k,10]
            rtime=source[k,12]
            #Where's the data
            strdepth='%.1f' % zs 
            syn_path=home+project_name+'/GFs/'+model_name+'_'+strdepth+'/'
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
            rtime=round(rtime/dt)*dt  #Rupture time has to be in the sampling rate of the GFs
            ess=tshift(ess,rtime)
            nss=tshift(nss,rtime)
            zss=tshift(zss,rtime)
            eds=tshift(eds,rtime)
            nds=tshift(nds,rtime)
            zds=tshift(zds,rtime)
            #get rake contribution and moment multiplier
            dsmult=sin(deg2rad(rake))
            ssmult=cos(deg2rad(rake))
            M=Mo/unitM
            if allclose(slip,0)==False:  #Only add things that matter
                print nfault+', SS='+str(ssmult)+', DS='+str(dsmult)+', Mscale='+str(M)
                #A'ight, add 'em up
                etotal=add_traces(ess,eds,ssmult,dsmult,M)
                ntotal=add_traces(nss,nds,ssmult,dsmult,M)
                ztotal=add_traces(zss,zds,ssmult,dsmult,M)
                #Add to previous subfault's results
                e=add_traces(e,etotal,1,1,1)
                n=add_traces(n,ntotal,1,1,1)
                z=add_traces(z,ztotal,1,1,1)
            else:
                print "No slip on subfault "+nfault+', ignoring it...'
                
        #Save results
        e[0].data=e[0].data.filled()  #This is a workaround of a bug in obspy
        n[0].data=n[0].data.filled()
        z[0].data=z[0].data.filled()
        e.write(outpath+sta+'.'+vord+'.e',format='SAC')
        n.write(outpath+sta+'.'+vord+'.n',format='SAC')
        z.write(outpath+sta+'.'+vord+'.z',format='SAC')
        

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
            dslength=source[k,10]
            #Where's the data
            strdepth='%.1f' % zs 
            syn_path=home+project_name+'/GFs/'+model_name+'_'+strdepth+'/'
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
    #Now extend both arrays and fill edges with 0's
    ss[0].trim(t1,t2,pad=True,fill_value=0)
    ds[0].trim(t1,t2,pad=True,fill_value=0)
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
    tstart=st[0].stats.starttime
    st[0].stats.starttime=tstart+tshift
    return st
        
        
    