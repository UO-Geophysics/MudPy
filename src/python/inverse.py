'''
Diego Melgar, 03.2014

Routines for solving dislocation ivnerse problems, static and dynamic.
'''
def getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,rupture_speeds,coord_type):
    '''
    Either load G from file or parse gflist file and assemble it from previous computations
    '''
    from numpy import genfromtxt,where,loadtxt,array,r_,concatenate,save,load
    from os import remove
    from os.path import split
    
    G_name=home+project_name+'/GFs/matrices/'+G_name
    if G_from_file==1: #load from file
        if G_name[-3:]!='npy':
            G_name=G_name+'.npy'
        print 'Loading G from file '+G_name
        G=load(G_name)
    else: #assemble G one data type at a time
        print 'Assembling G from synthetic computations...'
        #Read in GFlist and decide what to compute
        gf_file=home+project_name+'/data/station_info/'+GF_list
        mini_station=home+project_name+'/data/station_info/temp.sta'
        stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='S6')
        GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
        GFfiles=genfromtxt(gf_file,usecols=[8,9,10],dtype='S')
        #static field GFs
        kgf=2
        Gstatic=array([])
        if GF[:,kgf].sum()>0:
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,0])
                gftype='static'
                tdelay=0
                Gstatic=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay)
                remove(mini_station) #Cleanup  
        #Dispalcement waveform GFs
        kgf=3
        Gdisp=array([])
        if GF[:,kgf].sum()>0:
            #Load fault file
            source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,1])
                gftype='disp'
                for krup in range(len(rupture_speeds)):
                    tdelay=epi2subfault(epicenter,source,rupture_speeds[krup])
                    Gdisp_temp=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay)
                    if krup==0: #First rupture speed
                        Gdisp=Gdisp_temp
                    else:
                        Gdisp=r_[Gdisp,Gdisp_temp]
                remove(mini_station) #Cleanup 
        #Velocity waveforms
        kgf=4
        Gvel=array([])
        if GF[:,kgf].sum()>0:
            pass
        #Tsunami waveforms
        kgf=5
        Gtsun=array([])
        if GF[:,kgf].sum()>0:
            pass
        #Strain offsets
        kgf=6
        Gstrain=array([])
        if GF[:,kgf].sum()>0:
            pass
        #Done, now concatenate them all and save to pickle file
        G=concatenate([g for g in [Gstatic,Gdisp,Gvel,Gtsun,Gstrain] if g.size > 0])
        print 'Saving to '+G_name+' this might take just a second...'
        save(G_name,G)
    return G
    
    
def getdata(home,project_name,GF_list,rupture_speeds):
    '''
    Assemble the data vector
    '''
    from numpy import genfromtxt,where,array,append,r_,concatenate
    from obspy import read
    from forward import round_time
    

    #Read gf file and decide what needs tog et loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10],dtype='S')  
    stations=genfromtxt(gf_file,usecols=0,dtype='S')  
    #Read one column at a time
    kgf=0 #Static field
    dstatic=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
            dtemp=genfromtxt(GFfiles[i[ksta],kgf])
            e=dtemp[0]
            n=dtemp[0]
            u=dtemp[0]
            dstatic=append(dstatic,r_[n,e,u])
    #Displacements
    kgf=1
    ddisp=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        print 'Assembling displacement waveforms from '+stations[i[ksta]]+' into inversion.'
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        #Make sure they are rounded to a dt interval
        dt=e[0].stats.delta
        e[0].stats.starttime=round_time(e[0].stats.starttime,dt)
        n[0].stats.starttime=round_time(n[0].stats.starttime,dt)
        u[0].stats.starttime=round_time(u[0].stats.starttime,dt)
        ddisp=append(ddisp,r_[n[0].data,e[0].data,u[0].data])
    #If there is more than one rupture speed
    for k in range(1,len(rupture_speeds)):
        ddisp=r_[ddisp,ddisp]
    #Velocities
    kgf=2
    dvel=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        print 'Assembling displacement waveforms from '+stations[i[ksta]]+' into inversion.'
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        #Make sure they are rounded to a dt interval
        dt=e[0].stats.delta
        e[0].stats.starttime=round_time(e[0].stats.starttime,dt)
        n[0].stats.starttime=round_time(n[0].stats.starttime,dt)
        u[0].stats.starttime=round_time(u[0].stats.starttime,dt)
        dvel=append(dvel,r_[n[0].data,e[0].data,u[0].data])
    #If there is more than one rupture speed
    for k in range(1,len(rupture_speeds)):
        dvel=r_[dvel,dvel]
    #Tsunami
    kgf=3
    dtsun=array([])
    #Strain
    kgf=4
    dstrain=array([])            
    #Done, concatenate all and exit
    d=concatenate([dx for dx in [dstatic,ddisp,dvel,dtsun,dstrain] if dx.size > 0])
    return d
    
    
    
def getL(home,project_name,fault_name,bounds,regularization_type,nfaults):
    '''
    Make regularization matrix
    '''
    
    from numpy import loadtxt,zeros,ones
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    N=len(source) #No. of subfaults
    nstrike=nfaults[0]
    ndip=nfaults[1]
    #Initalize
    L=zeros((2*N,2*N))
    #Which L am I building?
    if regularization_type.lower()=='laplace':
        print 'Making discrete Laplace operator regularization matrix...'
        for kfault in range(N): #Loop over faults and fill regularization matrix
            stencil,correction=laplace_stencil(kfault,nstrike,ndip,bounds)
            #Add strike slip branches of stencil
            L[2*kfault,2*stencil]=1
            #Add dip slip branches of stencil
            L[2*kfault+1,2*stencil+1]=1
            #Add strike slip central node with correction
            L[2*kfault,2*kfault]=-4+correction
            #Add dip slip central node with correction
            L[2*kfault+1,2*kfault+1]=-4+correction
        return L
    else:
        print 'ERROR: Unknown regularization type ('+regularization_type+') requested.'
        return False
        
        
        
    
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

def laplace_stencil(ifault,nstrike,ndip,bounds):
    '''
    Find the index of the subfaults that make the laplacian stencil of fault number ifault
    
    Usage:
        stencil=laplace_stencil(ifault,nstrike,ndip)
        
    IN:
        ifault: subfault index number
        nstrike: number of along-strike fault segments
        ndip: number of along dip subfaults
        bounds: tuple containing information on the boundary conditions to be applied to the stencil
    OUT:
        stencil: indices of the subfaults that contribute to the laplacian
        correction: number to be ADDED to the central node to account for the point being on a corner/edge
    '''
    
    from numpy import array
    
    #What are the boundary condiitons
    top=bounds[0]
    bottom=bounds[1]
    left=bounds[2]
    right=bounds[3]
    row=ifault/nstrike #Row number corresponding to this subfault
    column=ifault-(nstrike*row)
    if row==0 and column==0: #Top right corner
        stencil=array([ifault+1,ifault+nstrike])
        if top.lower()=='free' and right.lower()=='free': #Both are free
            correction=2
        elif top.lower()=='free' or right.lower()=='free': #Only one is free
            correction=1
        else: #Both are locked
            correction=0
        return stencil,correction
    if row==0 and column==(nstrike-1): #Top left corner
        stencil=array([ifault-1,ifault+nstrike])
        if top.lower()=='free' and left.lower()=='free': #Both are free
            correction=2
        elif top.lower()=='free' or left.lower()=='free': #Only one is free
            correction=1
        else: #Both are locked
            correction=0
        return stencil,correction
    if row==(ndip-1) and column==0: #Bottom right corner
        stencil=array([ifault-nstrike,ifault+1])
        if bottom.lower()=='free' and right.lower()=='free': #Both are free
            correction=2
        elif bottom.lower()=='free' or right.lower()=='free': #Only one is free
            correction=1
        else: #Both are locked
            correction=0
        return stencil,correction
    if row==(ndip-1) and column==(nstrike-1): #Bottom left corner
        stencil=array([ifault-nstrike,ifault-1])
        if bottom.lower()=='free' and left.lower()=='free': #Both are free
            correction=2
        elif bottom.lower()=='free' or left.lower()=='free': #Only one is free
            correction=1
        else: #Both are locked
            correction=0
        return stencil,correction
    if row==0: #Top edge, NOT the corner
        stencil=array([ifault+1,ifault-1,ifault+nstrike])
        if top.lower()=='free': #Free boundary condition
            correction=1
        else: #Boundary is locked
            correction=0
        return stencil,correction
    if row==(ndip-1): #Bottom edge, NOT the corner
        stencil=array([ifault-1,ifault+1,ifault-nstrike])
        if bottom.lower()=='free': #Free boundary condition
            correction=1
        else: #Boundary is locked
            correction=0
        return stencil,correction
    if column==0: #Right edge, NOT the corner
        stencil=array([ifault-nstrike,ifault+nstrike,ifault+1])
        if right.lower()=='free': #Free boundary condition
            correction=1
        else: #Boundary is locked
            correction=0
        return stencil,correction
    if column==(nstrike-1): #left edge, NOT the corner
        stencil=array([ifault-nstrike,ifault+nstrike,ifault-1])
        if left.lower()=='free': #Free boundary condition
            correction=1
        else: #Boundary is locked
            correction=0
        return stencil,correction
    else: #Somewhere in the middle
        stencil=array([ifault-1,ifault+1,ifault-nstrike,ifault+nstrike])
        correction=0
        return stencil,correction



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
        
def mini_station_file(outfile,sta,lon,lat,gffiles):
    '''
    Make a temporary station file from a larger file
    '''
    f=open(outfile,'a')
    for k in range(len(sta)):
        out=sta[k]+'\t'+repr(round(lon[k],6))+'\t'+repr(round(lat[k],6))+'\t'+gffiles[k]+'\n'
        f.write(out)
    f.close()

    
            
def epi2subfault(epicenter,source,vr):
    '''
    Compute time delays from epicetner to subfault based on a give rupture speed
    
    Coordinates in lat/lon,depth(km). vr in km/s, returns tdelay in s
    '''
    from numpy import tile,sin,cos,deg2rad,sqrt
    #Compute distances from epi to subfault by converting to cartesian
    R=6371
    epicenter=tile(epicenter,(len(source),1))
    xepi=(R-epicenter[:,2])*sin(deg2rad(90-epicenter[:,1]))*cos(deg2rad(epicenter[:,0]))
    yepi=(R-epicenter[:,2])*sin(deg2rad(90-epicenter[:,1]))*sin(deg2rad(epicenter[:,0]))
    zepi=(R-epicenter[:,2])*cos(deg2rad(90-epicenter[:,1]))
    x=(R-source[:,3])*sin(deg2rad(90-source[:,2]))*cos(deg2rad(source[:,1]))
    y=(R-source[:,3])*sin(deg2rad(90-source[:,2]))*sin(deg2rad(source[:,1]))
    z=(R-source[:,3])*cos(deg2rad(90-source[:,2]))
    d=sqrt((xepi-x)**2+(yepi-y)**2+(zepi-z)**2)
    #Compute time associated with a given rupture speed
    tdelay=d/vr
    return tdelay   
    
    