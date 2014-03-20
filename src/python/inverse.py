'''
Diego Melgar, 03.2014

Routines for solving dislocation ivnerse problems, static and dynamic.
'''

def makeG(home,project_name,fault_name,model_name,station_file,gftype,tdelay,max_time):
    '''
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        station_file: File with coordinates of stations and data types
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
       
    OUT:
        Nothing
    '''
    from numpy import genfromtxt,loadtxt,zeros
    from string import rjust
    from obspy import read
    from forward import tshift,round_time
    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    Nsta=len(staname)
    #Initalize output variable
    G=zeros([len(staname)*3,Nfaults*2])
    #Loop over stations
    for ksta in range(Nsta):
            if gftype.lower()=='static': #Make matrix of static GFs
                #Where's the data
                syn_path=home+project_name+'/GFs/static/'
                #Loop over subfaults
                for kfault in range(Nfaults):
                    nfault='subfault'+rjust(str(int(source[kfault,0])),4,'0')
                    print 'Assembling static GFs for station '+staname[ksta]+' '+nfault
                    coseis_ss=loadtxt(syn_path+staname[ksta]+'.'+nfault+'.SS.static.enu')
                    ess=coseis_ss[0]
                    nss=coseis_ss[1]
                    zss=coseis_ss[2]
                    coseis_ds=loadtxt(syn_path+staname[ksta]+'.'+nfault+'.DS.static.enu')
                    eds=coseis_ds[0]
                    nds=coseis_ds[1]
                    zds=coseis_ds[2]
                    #Place into G matrix
                    G[3*ksta,2*kfault]=nss   ; G[3*ksta,2*kfault+1]=nds    #North
                    G[3*ksta+1,2*kfault]=ess ; G[3*ksta+1,2*kfault+1]=eds  #East
                    G[3*ksta+2,2*kfault]=zss ; G[3*ksta+2,2*kfault+1]=zds  #Up
            if gftype.lower()=='disp' or gftype.lower=='vel':  #Full waveforms
                if gftype.lower()=='disp':
                    vord='disp'
                else:
                    vord='vel'
                #determine length of waveforms
                #Loop over subfaults
                for kfault in range(Nfaults):
                    print 'Assembling '+vord+' waveform GFs for station '+staname[ksta]+' '+nfault
                    #Get subfault GF directory
                    nsub='sub'+rjust(str(int(source[kfault,0])),4,'0')
                    nfault='subfault'+rjust(str(int(source[kfault,0])),4,'0')
                    strdepth='%.4f' % source[kfault,3]
                    syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                    #Get synthetics
                    ess=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
                    nss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
                    zss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
                    eds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
                    nds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
                    zds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
                    #Time shift them according to subfault rupture time, zero pad, round to dt interval
                    #and extend to maximum time
                    dt=ess[0].stats.delta
                    ess=tshift(ess,tdelay[kfault])
                    ess[0].stats.starttime=round_time(ess[0].stats.starttime,dt)
                    ess=pad2zero(ess)
                    ess=padend(ess,max_time,vord)
                    nss=tshift(nss,tdelay[kfault])
                    nss[0].stats.starttime=round_time(nss[0].stats.starttime,dt)
                    nss=pad2zero(nss)
                    nss=padend(nss,max_time,vord)
                    zss=tshift(zss,tdelay[kfault])
                    zss[0].stats.starttime=round_time(zss[0].stats.starttime,dt)
                    zss=pad2zero(zss)
                    zss=padend(zss,max_time,vord)
                    eds=tshift(eds,tdelay[kfault])
                    eds[0].stats.starttime=round_time(eds[0].stats.starttime,dt)
                    eds=pad2zero(eds)
                    eds=padend(eds,max_time,vord)
                    nds=tshift(nds,tdelay[kfault])
                    nds[0].stats.starttime=round_time(nds[0].stats.starttime,dt)
                    nds=pad2zero(nds)
                    nds=padend(nds,max_time,vord)
                    zds=tshift(zds,tdelay[kfault])
                    zds[0].stats.starttime=round_time(zds[0].stats.starttime,dt)
                    zds=pad2zero(zds)
                    zds=padend(zds,max_time,vord)
                    #Initalize G
                    if kfault==0 and ksta==0: 
                        npts=ess[0].stats.npts
                        G=zeros([3*Nsta*npts,2*Nfaults])
                    #Insert them into G
                    G[3*ksta*npts:3*ksta*npts+npts,2*kfault]=nss[0].data
                    G[3*ksta*npts:3*ksta*npts+npts,2*kfault+1]=nds[0].data
                    G[3*ksta*npts+npts:3*ksta*npts+2*npts,2*kfault]=ess[0].data
                    G[3*ksta*npts+npts:3*ksta*npts+2*npts,2*kfault+1]=eds[0].data
                    G[3*ksta*npts+2*npts:3*ksta*npts+3*npts,2*kfault]=zss[0].data
                    G[3*ksta*npts+2*npts:3*ksta*npts+3*npts,2*kfault+1]=zds[0].data 
            if gftype.lower()=='tsun':
                pass                
            if gftype.lower()=='strain':
                pass                                
    return G
    
    
#==================              Random Tools            ======================
  
def pad2zero(st):
    '''
    Extend a stream object to zero time and pad with 0's
    
    Usage:
        st=pad2zero(st)
    ''' 
    from numpy import zeros,r_
    from obspy import Stream,Trace
    
    #Figure out how many zeros I will need
    npad=st[0].stats.starttime.minute*60+st[0].stats.starttime.second+st[0].stats.starttime.microsecond/1e6
    npad=npad/st[0].stats.delta
    #Build zeros trace
    z=Stream(Trace())
    z[0].stats.starttime=st[0].stats.starttime
    z[0].stats.starttime.minute=0
    z[0].stats.starttime.second=0
    z[0].stats.starttime.microsecond=0
    z[0].stats.delta=st[0].stats.delta
    z[0].stats.npts=npad
    z[0].data=zeros(npad)
    #Add and make output stream
    st_out=Stream(Trace())
    st_out[0].stats.starttime=z[0].stats.starttime
    st_out[0].stats.delta=st[0].stats.delta
    st_out[0].stats.npts=z[0].stats.npts+st[0].stats.npts
    st_out[0].data=r_[z[0].data,st[0].data]
    
    return st_out
    
def padend(st,tend,vord):
    '''
    Extend a waveform to alonger time and fill with a given value
    '''
    st_out=st.copy()
    if vord.lower()=='disp': #Fillw ith mean of last 100 samples
        fillend=st[0].data[-100:-1].mean()
    else: #Fill with zeros
        fillend=0
    st_out[0].trim(st[0].stats.starttime,tend,pad=True,fill_value=fillend)
    return st_out