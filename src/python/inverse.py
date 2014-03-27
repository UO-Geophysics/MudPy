'''
Diego Melgar, 03.2014

Routines for solving dislocation ivnerse problems, static and dynamic.
'''

def makeG(home,project_name,fault_name,model_name,station_file,gftype,tdelay):
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
    from numpy import genfromtxt,loadtxt,zeros,r_
    from string import rjust
    from obspy import read
    from forward import tshift,round_time
    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    datafiles=genfromtxt(station_file,dtype="S",usecols=3)
    Nsta=len(staname)
    insert_position=0
    #Loop over stations
    for ksta in range(Nsta):
            if gftype.lower()=='static': #Make matrix of static GFs
                #Initalize output variable
                Gtemp=zeros([3,Nfaults*2])
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
                    Gtemp[0,2*kfault]=nss   ; Gtemp[0,2*kfault+1]=nds    #North
                    Gtemp[1,2*kfault]=ess ; Gtemp[1,2*kfault+1]=eds  #East
                    Gtemp[2,2*kfault]=zss ; Gtemp[2,2*kfault+1]=zds  #Up
                    #Append to G
                if ksta==0: #First station, create array 
                    G=Gtemp
                else: #Just append
                    G=r_[G,Gtemp]
            if gftype.lower()=='disp' or gftype.lower=='vel':  #Full waveforms
                if gftype.lower()=='disp':
                    vord='disp'
                else:
                    vord='vel'
                #Loop over subfaults
                for kfault in range(Nfaults):
                    #Get subfault GF directory
                    nsub='sub'+rjust(str(int(source[kfault,0])),4,'0')
                    nfault='subfault'+rjust(str(int(source[kfault,0])),4,'0')
                    strdepth='%.4f' % source[kfault,3]
                    syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                    print 'Assembling '+vord+' waveform GFs for station '+staname[ksta]+' '+nfault
                    #Get synthetics
                    ess=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
                    nss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
                    zss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
                    eds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
                    nds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
                    zds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
                    #Load raw data, this will be used to trim the GFs
                    edata=read(datafiles[ksta]+'.e')
                    ndata=read(datafiles[ksta]+'.n')
                    udata=read(datafiles[ksta]+'.u')
                    #Time shift them according to subfault rupture time, zero pad, round to dt interval
                    #and extend to maximum time
                    dt=ess[0].stats.delta
                    ess=tshift(ess,tdelay[kfault])
                    ess[0].stats.starttime=round_time(ess[0].stats.starttime,dt)
                    ess=prep_synth(ess,edata)
                    nss=tshift(nss,tdelay[kfault])
                    nss[0].stats.starttime=round_time(nss[0].stats.starttime,dt)
                    nss=prep_synth(nss,ndata)
                    zss=tshift(zss,tdelay[kfault])
                    zss[0].stats.starttime=round_time(zss[0].stats.starttime,dt)
                    zss=prep_synth(zss,udata)
                    eds=tshift(eds,tdelay[kfault])
                    eds[0].stats.starttime=round_time(eds[0].stats.starttime,dt)
                    eds=prep_synth(eds,edata)
                    nds=tshift(nds,tdelay[kfault])
                    nds[0].stats.starttime=round_time(nds[0].stats.starttime,dt)
                    nds=prep_synth(nds,ndata)
                    zds=tshift(zds,tdelay[kfault])
                    zds[0].stats.starttime=round_time(zds[0].stats.starttime,dt)
                    zds=prep_synth(zds,udata)
                    #Insert into Gtemp then append to G
                    if kfault==0 and ksta==0: #It's the first subfault and station, initalize G
                        G=gdims(datafiles,Nfaults) #Survey all stations to decide size of G
                    if kfault==0: #Initalize Gtemp
                        npts=edata[0].stats.npts
                        Gtemp=zeros([3*npts,Nfaults*2])      
                    #Insert synthetics into Gtemp
                    Gtemp[0:npts,2*kfault]=nss[0].data
                    Gtemp[0:npts,2*kfault+1]=nds[0].data
                    Gtemp[npts:2*npts,2*kfault]=ess[0].data
                    Gtemp[npts:2*npts,2*kfault+1]=eds[0].data
                    Gtemp[2*npts:3*npts,2*kfault]=zss[0].data
                    Gtemp[2*npts:3*npts,2*kfault+1]=zds[0].data
                #After looping through all subfaults Insert Gtemp into G
                G[insert_position:insert_position+3*npts,:]=Gtemp
                insert_position+=3*npts #Update for next station
            if gftype.lower()=='tsun':
                pass                
            if gftype.lower()=='strain':
                pass                                
    return G
    
    
#==================              Random Tools            ======================
      
def prep_synth(syn,st):
    '''
    Extend syntetic to start time of data and cut it to end time of data, make sure
    synthetic ALWAYS ends after data
    '''
    from numpy import zeros,r_
    #What's the difference ins tart times?
    t1syn=syn[0].stats.starttime
    t1st=st[0].stats.starttime
    dt=t1syn-t1st
    if dt>=0: #Synthetic starts before data, pad with zeros
        #How many zeros do I need
        npad=dt/st[0].stats.delta
        z=zeros(npad)
        syn[0].data=r_[z,syn[0].data]
        syn[0].stats.starttime=t1st
    else: #Synthetic starts after the waveform, crop to start fo waveform
        syn[0].trim(t1st)
    #Now deal with end times
    t2syn=syn[0].stats.endtime
    t2st=st[0].stats.endtime
    dt=t2syn-t2st
    if dt>=0: #Synthetic ends after data, crop it
        syn[0].trim(endtime=t2st)
    else: #Syntetic ends before data, throw an error
        print "ERROR: Synthetic end time is before data end time, recompute longer syntehtics please."
        return 'Error in GF length'
    return syn
        
def gdims(datafiles,nfaults):
    '''
    Survey the data files to determine what dimension G will be and return a matrix of zeros 
    with the required dimensions
    '''
    
    from obspy import read
    from numpy import zeros
    
    npts=0
    for k in range(len(datafiles)):
        e=read(datafiles[k]+'.e')
        n=read(datafiles[k]+'.n')
        u=read(datafiles[k]+'.u')
        if e[0].stats.npts==n[0].stats.npts==u[0].stats.npts:
            npts+=e[0].stats.npts
        else:
            print str(e[0].stats.npts)+' pts in east component'
            print str(n[0].stats.npts)+' pts in north component'
            print str(u[0].stats.npts)+' pts in up component'
            print 'ERROR: The 3 components of data are not the same length'
            return 'Error in forming G'
    G=zeros([3*npts,nfaults*2])
    return G            
        
    
    
    