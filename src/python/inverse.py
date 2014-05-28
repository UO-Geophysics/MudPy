'''
Diego Melgar, 03.2014

Routines for solving dislocation inverse problems, static and kinematic.

Functions in this module:
    
    * getG() - Assembles Green functions matrix for ALL data types
    * makeG() - Assembles Green fucntions for a particular data type
    * getdata() - Assembles data vector for inversion
    * getLs() - Assembles spatial regularization matrix based on finite difference Laplacian
    * getLt() - Assembles temporal regularization using finite difference first derivatives
    * get_data_weights() - Assemble matrix of data weights

'''



def getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,rupture_speed,
        num_windows,coord_type,decimate):
    '''
    Assemble Green functions matrix. If requested will parse all available synthetics on file and build the matrix.
    Otherwise, if it exists, it will be loaded from file 
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        GF_list: Name of GF control file
        G_from_file: if =0 build G from synthetics on file. If =1 then load from file
        G_name: If building G fromsynthetics then this is the name G will be saved with
            in binary .npy format. If loading from file this is the name to be looked for. 
            It is not necessary to supply the .npy extension
        epicenter: Epicenter coordinates
        rupture_speed: Fastest rupture speed allowed in the problem
        num_windows: Number of temporal rupture windows allowed
        coord_type: =0 for cartesian, =1 for lat/lon
        decimate: Constant decimationf actor applied to GFs, set =0 for no decimation
        
    OUT:
        G: Fully assembled GF matrix
    '''
    
    from numpy import arange,genfromtxt,where,loadtxt,array,c_,concatenate,save,load,dot
    from os import remove
    from os.path import split
    
    G_name=home+project_name+'/GFs/matrices/'+G_name
    K_name=G_name+'.K'
    if G_from_file==1: #load from file
        if G_name[-3:]!='npy':
            K_name=K_name+'.npy'
            G_name=G_name+'.npy'
        print 'Loading G from file '+G_name
        G=load(G_name)
        K=load(K_name)
    else: #assemble G one data type at a time
        print 'Assembling G from synthetic computations...'
        #Read in GFlist and decide what to compute
        gf_file=home+project_name+'/data/station_info/'+GF_list
        mini_station=home+project_name+'/data/station_info/tempG.sta'
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
                Gstatic=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay,decimate)
                remove(mini_station) #Cleanup  
        #Dispalcement waveform GFs
        kgf=3
        Gdisp=array([])
        if GF[:,kgf].sum()>0:
            #Load fault file
            source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
            trise=source[0,7]
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,1])
                gftype='disp'
                #Decide on delays for each time window (50% overlap)
                trupt=arange(0,num_windows)*trise/2
                for krup in range(num_windows):
                    print 'Working on window '+str(krup+1)
                    tdelay=epi2subfault(epicenter,source,rupture_speed,trupt[krup])
                    Gdisp_temp=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay,decimate)
                    if krup==0: #First rupture speed
                        Gdisp=Gdisp_temp
                    else:
                        Gdisp=c_[Gdisp,Gdisp_temp]
                remove(mini_station) #Cleanup 
        #Velocity waveforms
        kgf=4
        Gvel=array([])
        if GF[:,kgf].sum()>0:
            #Load fault file
            source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
            trise=source[0,7]
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,1])
                gftype='vel'
                #Decide on delays for each time window (50% overlap)
                trupt=arange(0,num_windows)*trise/2
                for krup in range(num_windows):
                    print 'Working on window '+str(krup+1)
                    tdelay=epi2subfault(epicenter,source,rupture_speed,trupt[krup])
                    Gvel_temp=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay,decimate)
                    if krup==0: #First rupture speed
                        Gvel=Gvel_temp
                    else:
                        Gvel=c_[Gvel,Gvel_temp]
                remove(mini_station) #Cleanup 
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
        #Done, now concatenate them all ccompute transpose product and save
        G=concatenate([g for g in [Gstatic,Gdisp,Gvel,Gtsun,Gstrain] if g.size > 0])
        #K=dot(G.T,G)
        print 'Computing G\'G and saving to '+G_name+' this might take just a second...'
        save(G_name,G)
        #save(K_name,K)
    return G
    
def makeG(home,project_name,fault_name,model_name,station_file,gftype,tdelay,decimate):
    '''
    This routine is called from getG and will assemble the GFs from available synthetics
    depending on data type requested (statics, dispalcement or velocity waveforms).
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        station_file: File with coordinates of stations and data types
        gftype: ='static' if assembling static field GFs, ='disp' if assembling displacement
            waveforms. ='vel' if assembling velocity waveforms.
        tdelay: Vector of delay times to be applied to each time window
        decimate: Constant decimationf actor applied to GFs, set =0 for no decimation
       
    OUT:
        G: Partially assembled GF with all synthetics from a particular data type
    '''
    from numpy import genfromtxt,loadtxt,zeros
    from string import rjust
    from obspy import read
    from forward import tshift
    from green import stdecimate
    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    datafiles=genfromtxt(station_file,dtype="S",usecols=3)
    Nsta=len(staname)
    insert_position=0
    #Initalize G for faster assignments
    if gftype.lower()=='static': #Initialize output matrix
        G=zeros((Nsta*3,Nfaults*2))
    else:
        pass #For disp or vel waveforms G is initalized below
    #Loop over stations
    for ksta in range(Nsta):
            if gftype.lower()=='static': #Make matrix of static GFs
                print 'Assembling static GFs for station '+staname[ksta]
                #Initalize output variable
                Gtemp=zeros([3,Nfaults*2])
                #Where's the data
                syn_path=home+project_name+'/GFs/static/'
                #Loop over subfaults
                for kfault in range(Nfaults):
                    if kfault%10==0:
                        print '... working on subfault '+str(kfault)+' of '+str(Nfaults)
                    nfault='subfault'+rjust(str(int(source[kfault,0])),4,'0')
                    coseis_ss=loadtxt(syn_path+staname[ksta]+'.'+nfault+'.SS.static.neu')
                    nss=coseis_ss[0]
                    ess=coseis_ss[1]
                    zss=coseis_ss[2]
                    coseis_ds=loadtxt(syn_path+staname[ksta]+'.'+nfault+'.DS.static.neu')
                    nds=coseis_ds[0]
                    eds=coseis_ds[1]
                    zds=coseis_ds[2]
                    #Place into G matrix
                    Gtemp[0,2*kfault]=nss   ; Gtemp[0,2*kfault+1]=nds    #North
                    Gtemp[1,2*kfault]=ess ; Gtemp[1,2*kfault+1]=eds  #East
                    Gtemp[2,2*kfault]=zss ; Gtemp[2,2*kfault+1]=zds  #Up
                    #Append to G
                #Append to output matrix
                G[ksta*3:ksta*3+3,:]=Gtemp
            if gftype.lower()=='disp' or gftype.lower()=='vel':  #Full waveforms
                if gftype.lower()=='disp':
                    vord='disp'
                else:
                    vord='vel'
                print 'Assembling '+vord+' waveform GFs for station '+staname[ksta]
                #Loop over subfaults
                for kfault in range(Nfaults):
                    if kfault%10==0:
                        print '... working on subfault '+str(kfault)+' of '+str(Nfaults)
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
                    #Time shift them according to subfault rupture time, zero pad, round to dt interval,decimate
                    #and extend to maximum time
                    ess=tshift(ess,tdelay[kfault])
                    nss=tshift(nss,tdelay[kfault])
                    zss=tshift(zss,tdelay[kfault])
                    eds=tshift(eds,tdelay[kfault])
                    nds=tshift(nds,tdelay[kfault])
                    zds=tshift(zds,tdelay[kfault])
                    if decimate>0: 
                        ess=stdecimate(ess,decimate)
                        nss=stdecimate(nss,decimate)
                        zss=stdecimate(zss,decimate)    
                        eds=stdecimate(eds,decimate)
                        nds=stdecimate(nds,decimate)
                        zds=stdecimate(zds,decimate)
                    if decimate>0 and kfault==0:#Only need to do this once
                        #Load raw data, this will be used to trim the GFs
                        edata=read(datafiles[ksta]+'.e')
                        ndata=read(datafiles[ksta]+'.n')
                        udata=read(datafiles[ksta]+'.u')
                        edata=stdecimate(edata,decimate)
                        ndata=stdecimate(ndata,decimate)
                        udata=stdecimate(udata,decimate)
                        #How many points left in the tiem series
                        npts=edata[0].stats.npts
                        print "... "+str(npts)+" data points left over after decimation"
                    #Now time align stuff (This is where we might have some timing issues, check later, consider re-sampling to data time vector)                                       
                    ess=resample_to_data(ess,edata)
                    ess=prep_synth(ess,edata)
                    nss=resample_to_data(nss,ndata)
                    nss=prep_synth(nss,ndata)
                    zss=resample_to_data(zss,udata)
                    zss=prep_synth(zss,udata)
                    eds=resample_to_data(eds,edata)
                    eds=prep_synth(eds,edata)
                    nds=resample_to_data(nds,ndata)
                    nds=prep_synth(nds,ndata)
                    zds=resample_to_data(zds,udata)
                    zds=prep_synth(zds,udata)
                    #Insert into Gtemp then append to G
                    if kfault==0 and ksta==0: #It's the first subfault and station, initalize G
                        G=gdims(datafiles,Nfaults,decimate) #Survey all stations to decide size of G
                    if kfault==0: #Initalize Gtemp (different size for each station)
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
      
      
def getdata(home,project_name,GF_list,decimate):
    '''
    Assemble the data vector for all data types
    
    IN:
        home: Home directory
        project_name: Name of the problem
        GF_list: Name of GF control file
        decimate: Constant decimationf actor applied to GFs, set =0 for no decimation
        
    OUT:
        d: The data vector
    '''
    from numpy import genfromtxt,where,array,append,r_,concatenate,zeros
    from obspy import read
    from forward import round_time
    from green import stdecimate
    

    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10],dtype='S')  
    stations=genfromtxt(gf_file,usecols=0,dtype='S')  
    #Read one column at a time
    kgf=0 #Static field
    dstatic=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
            print 'Assembling static offsets from '+stations[i[ksta]]+' into data vector.'
            dtemp=genfromtxt(GFfiles[i[ksta],kgf])
            n=dtemp[0]
            e=dtemp[1]
            u=dtemp[2]
            dstatic=append(dstatic,r_[n,e,u])
    #Displacements
    kgf=1
    ddisp=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        print 'Assembling displacement waveforms from '+stations[i[ksta]]+' into data vector.'
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        #Decimate
        if decimate>0:
            n=stdecimate(n,decimate)
            e=stdecimate(e,decimate)
            u=stdecimate(u,decimate)
        #Make sure they are rounded to a dt interval
        dt=e[0].stats.delta
        e[0].stats.starttime=round_time(e[0].stats.starttime,dt)
        n[0].stats.starttime=round_time(n[0].stats.starttime,dt)
        u[0].stats.starttime=round_time(u[0].stats.starttime,dt)
        ddisp=append(ddisp,r_[n[0].data,e[0].data,u[0].data])
    #Velocities
    kgf=2
    dvel=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        print 'Assembling velocity waveforms from '+stations[i[ksta]]+' into data vector.'
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        #Decimate
        if decimate>0:
            n=stdecimate(n,decimate)
            e=stdecimate(e,decimate)
            u=stdecimate(u,decimate)
        #Make sure they are rounded to a dt interval
        dt=e[0].stats.delta
        e[0].stats.starttime=round_time(e[0].stats.starttime,dt)
        n[0].stats.starttime=round_time(n[0].stats.starttime,dt)
        u[0].stats.starttime=round_time(u[0].stats.starttime,dt)
        dvel=append(dvel,r_[n[0].data,e[0].data,u[0].data])
    #Tsunami
    kgf=3
    dtsun=array([])
    #Strain
    kgf=4
    dstrain=array([])            
    #Done, concatenate all, convert to column vector and exit
    d=concatenate([dx for dx in [dstatic,ddisp,dvel,dtsun,dstrain] if dx.size > 0])
    D=zeros((d.shape[0],1))
    D[:,0]=d
    return D
          
    
def getLs(home,project_name,fault_name,nfaults,num_windows,bounds):
    '''
    Make spatial regularization matrix based on finite difference Lapalce operator.
    This routine will request adjustments depending on the boundary conditions requested
    on the edges of the fault model.
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        nfaults: Total number of faults in the model
        num_windows: Number of rupture windows
        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
            Possible values for each element of the tuple are 'free' for a free boundary condition
            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
            edge free would be bounds=('free', 'locked', 'locked', 'locked')

    OUT:
        Lout: The regularization matrix
    '''
    
    from numpy import loadtxt,zeros,tile
    from scipy.linalg import block_diag
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    N=len(source) #No. of subfaults
    nstrike=nfaults[0]
    ndip=nfaults[1]
    #Initalize
    L=zeros((2*N,2*N))
    #Which L am I building?
    print 'Making discrete Laplace operator regularization matrix...'
    for kfault in range(N):#Loop over faults and fill regularization matrix
        stencil,values=laplace_stencil(kfault,nstrike,ndip,bounds)
        #Add strike slip branches of stencil
        L[2*kfault,2*stencil]=values
        #Add dip slip branches of stencil
        L[2*kfault+1,2*stencil+1]=values
    if num_windows==1: #Only one rupture speed
        Lout=L 
    else: #Multiple rupture speeds, smooth total moment laplacian
        Lout=L
        Lout=tile(Lout,(1,num_windows))/num_windows
        #for k in range(num_windows-1):
        #    Lout=block_diag(Lout,L)
    return Lout
        
        
def getLt(home,project_name,fault_name,num_windows):
    '''
    Make temporal regularization matrix using forward differences for windows 1 
    through N-1 and backwards differences for window N
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        num_windows: Number of temporal slip windows
        
    OUT:
        L: A square matrix of derivatives
    '''
    
    from numpy import loadtxt,zeros,eye,arange
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    N=len(source) #No. of subfaults
    if num_windows<2: #Duh
        print 'WARNING temporal regularization is unecessary when only employing 1 time window. Returning zeros.'
        Lout=zeros((2*N,2*N))
        return Lout
    #Initalize
    L=eye(2*N*num_windows)
    #Ltest=zeros((N*num_windows,2*N*num_windows))
    print 'Making first derivative temporal regularization matrix...'
    #Forward difference indices
    iforward=arange(N*2*(num_windows-1))
    L[iforward,iforward+(2*N)]=-1
    #
    #i1=arange(0,N*num_windows)
    #i2=arange(0,2*N*num_windows,2)
    #Ltest[i1,i2]=1
    #Ltest[i1,i2+1]=1
    #i1=arange(0,N*num_windows-N)
    #i2=arange(2*N,2*N*num_windows,2)
    #Ltest[i1,i2]=-1
    #Ltest[i1,i2+1]=-1
    #Ltest=Ltest/2
    #return Ltest
    #
    #Backwards differences for last window
    #iback=arange(N*2*(num_windows-1),N*2*num_windows)
    #L[iback,iback-(2*N)]=-1
    return L


def get_data_weights(home,project_name,GF_list,d,decimate):
    '''
    Assemble matrix of data weights from sigmas of observations
    '''    
    from numpy import genfromtxt,where,zeros,ones
    from obspy import read

    print 'Computing data weights...'
    #Read gf file and decide what needs tog et loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10],dtype='S')
    weights=genfromtxt(gf_file,usecols=range(13,28),dtype='f')
    #Initalize
    w=zeros(len(d))
    kinsert=0
    #Deal with decimation
    if decimate==0:
        decimate=1
    #Static weights
    kgf=0
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        w[kinsert]=1/weights[i[ksta],0] #North
        w[kinsert+1]=1/weights[i[ksta],1] #East
        w[kinsert+2]=1/weights[i[ksta],2] #Up
        kinsert=kinsert+3
    #Displacement waveform weights
    kgf=1
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        #Read waveform to determine length of insert
        st=read(GFfiles[i[ksta],kgf]+'.n')
        nsamples=st[0].stats.npts/decimate
        wn=(1/weights[i[ksta],3])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=wn
        kinsert=kinsert+nsamples
        we=(1/weights[i[ksta],4])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=we
        kinsert=kinsert+nsamples
        wu=(1/weights[i[ksta],5])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=wu
        kinsert=kinsert+nsamples
    #velocity waveform weights
    kgf=2
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        #Read waveform to determine length of insert
        st=read(GFfiles[i[ksta],kgf]+'.n')
        nsamples=st[0].stats.npts/decimate
        wn=(1/weights[i[ksta],6])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=wn
        kinsert=kinsert+nsamples
        we=(1/weights[i[ksta],7])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=we
        kinsert=kinsert+nsamples
        wu=(1/weights[i[ksta],8])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=wu
        kinsert=kinsert+nsamples
    #Tsunami
    kgf=3
    #Strain
    kgf=4
    #Make W and exit
    return w

    
    
#=================        Write inversion results      =========================
    
def write_model(home,project_name,run_name,fault_name,model_name,rupture_speed,num_windows,epicenter,sol,num):
    '''
    Write inversion results to .inv file
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        run_name: Name of inversion run
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        rupture_speed: Fastest rupture speed allowed in the problem
        num_windows: Number of temporal rupture windows allowed
        epicenter: Epicenter coordinates
        sol: The solution vector from the inversion
        num: ID number of the inversion
        GF_list: Name of GF control file
    OUT:
        Nothing
    '''
    
    from numpy import genfromtxt,loadtxt,arange,zeros,c_,savetxt,r_
    from forward import get_mu
    from string import rjust
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    trise=f[0,7]
    #Open structure file
    mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
    #Get slip quantities
    iss=2*arange(len(f)*num_windows)
    ids=2*arange(len(f)*num_windows)+1
    ss=sol[iss]
    ds=sol[ids]
    #Get rigidities
    mu=zeros(len(ds))
    trup=zeros(len(ds))
    j=0
    for krup in range(num_windows):
        for k in range(len(f)):
            mu[j]=get_mu(mod,f[k,3])
            j+=1
    #Get rupture start times
    trupt=arange(0,num_windows)*trise/2 #Time delays fore ach sub-window
    for krup in range(num_windows):
        trup[krup*(len(ds)/num_windows):(krup+1)*(len(ds)/num_windows)]=epi2subfault(epicenter,f,rupture_speed,trupt[krup])
    #Prepare for output
    out1=f[:,0:8]
    out2=f[:,8:10]
    for k in range(num_windows-1):
        out1=r_[out1,f[:,0:8]]
        out2=r_[out2,f[:,8:10]]
    out=c_[out1,ss,ds,out2,trup,mu]
    outdir=home+project_name+'/output/inverse_models/models/'+run_name+'.'+rjust(str(num),4,'0')+'.inv'
    #CHANGE this to rupture definition as #No  x            y        z(km)      str     dip      rake       rise    dura     slip    ss_len  ds_len rupt_time
    fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
    print '... writing model results to file '+outdir
    savetxt(outdir,out,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
        
    
    
def write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,num):
    '''
    Output synthetics as sac for displacement or velocity waveforms and ascii for static field
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        run_name: Name of inversion run
        sol: The solution vector from the inversion
        ds: The predicted data ds=G*m
        num: ID number of the inversion
    OUT:
        Nothing
    '''
    
    from obspy import read
    from numpy import array,savetxt,where,genfromtxt
    from string import rjust
    
    print '... computing and saving synthetics...'
    num=rjust(str(num),4,'0')
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='S')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10],skip_header=1,dtype='S')
    #Separate into its constituent parts (statics,displacaments, velocities, etc...)
    kinsert=0
    #Statics
    kgf=0
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            neu=array([ds[kinsert],ds[kinsert+1],ds[kinsert+2]])
            kinsert+=3
            savetxt(home+project_name+'/output/inverse_models/statics/'+run_name+'.'+num+'.'+sta+'.static.neu',neu)
    #Displacement
    kgf=1
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            n=read(GFfiles[i[ksta],kgf]+'.n')
            e=read(GFfiles[i[ksta],kgf]+'.e')
            u=read(GFfiles[i[ksta],kgf]+'.u')
            npts=n[0].stats.npts
            n[0].data=ds[kinsert:kinsert+npts]
            e[0].data=ds[kinsert+npts:kinsert+2*npts]
            u[0].data=ds[kinsert+2*npts:kinsert+3*npts]
            kinsert+=3*npts
            n.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.n.sac',format='SAC')
            e.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.e.sac',format='SAC')
            u.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.u.sac',format='SAC')
    #Velocity
    kgf=2
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            n=read(GFfiles[i[ksta],kgf]+'.n')
            e=read(GFfiles[i[ksta],kgf]+'.e')
            u=read(GFfiles[i[ksta],kgf]+'.u')
            npts=n[0].stats.npts
            n[0].data=ds[kinsert:kinsert+npts]
            e[0].data=ds[kinsert+npts:kinsert+2*npts]
            u[0].data=ds[kinsert+2*npts:kinsert+3*npts]
            kinsert+=3*npts
            n.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.n.sac',format='SAC')
            e.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.e.sac',format='SAC')
            u.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.u.sac',format='SAC')
            
        
def write_log(home,project_name,run_name,k,rupture_speed,num_windows,lambda_spatial,lambda_temporal,
        beta,L2,Lm,VR,ABIC,Mo,Mw,velmod,fault,g_name,gflist,solver):
    '''
    Write inversion sumamry to .log file
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        run_name: Name of inversion run
        k: Inversion run number
        rupture_speed: Fastest rupture speed allowed
        num_windows: Number of temporal rupture windows
        lambda_spatial: Spatial regularization parameter
        lambda_temporal: Temporal regularization parameter
        beta: Angular offset applied to rake
        L2: L2 norm of ivnersion L2=||Gm-d||
        Lm: Model norm Lm=||L*m||
        VR: Variance reduction
        ABIC: Value of Akaike's Bayesian ifnormation criterion
        Mo: Moment in N-m
        Mw: Moment magnitude
        velmod: Earth structure model used
        fault: Fault model used
        g_name: GF matrix used
        gflist: GF control file sued
        solver: Type of solver used
    OUT:
        Nothing
    '''
    
    from string import rjust
    
    num=rjust(str(k),4,'0')
    f=open(home+project_name+'/output/inverse_models/models/'+run_name+'.'+num+'.log','w')
    f.write('Project: '+project_name+'\n')
    f.write('Run name: '+run_name+'\n')
    f.write('Run number: '+num+'\n')
    f.write('Velocity model: '+velmod+'\n')
    f.write('Fault model: '+fault+'\n')
    f.write('G name: '+g_name+'\n')
    f.write('GF list: '+gflist+'\n')
    f.write('Solver: '+solver+'\n')
    f.write('lambda_spatial = '+repr(lambda_spatial)+'\n')
    f.write('lambda_temporal = '+repr(lambda_temporal)+'\n')
    f.write('Beta(degs) = '+repr(beta)+'\n')
    f.write('Mean rupture velocity (km/s) = '+str(rupture_speed)+'\n')
    f.write('Number of rupture windows = '+str(num_windows)+'\n')
    f.write('L2 = '+repr(L2)+'\n')
    f.write('VR(%) = '+repr(VR)+'\n')
    f.write('Lm = '+repr(Lm)+'\n')
    f.write('ABIC = '+repr(ABIC)+'\n')
    f.write('M0(N-m) = '+repr(Mo)+'\n')
    f.write('Mw = '+repr(Mw)+'\n')
    f.close()
    
    
    
#==================              Random Tools            ======================

def laplace_stencil(ifault,nstrike,ndip,bounds):
    '''
    Find the index of the subfaults that make the laplacian stencil of fault number ifault
    It assumes all boundaries are initally locked. After assigning stencil values it parses
    the variable 'bounds' and makes corrections if any boundaries were requested to be
    'free'.
    
    Usage:
        stencil=laplace_stencil(ifault,nstrike,ndip)
        
    IN:
        ifault: subfault index number
        nstrike: number of along-strike fault segments
        ndip: number of along dip subfaults
        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
            Possible values for each element of the tuple are 'free' for a free boundary condition
            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
            edge free would be bounds=('free', 'locked', 'locked', 'locked')
    OUT:
        stencil: indices of the subfaults that contribute to the laplacian
        values: nuemrical values of the stencil
    '''
    
    from numpy import array
    
    #Get boundary conditions
    top=bounds[0]
    bottom=bounds[1]
    left=bounds[2]
    right=bounds[3]
    #Create stencil
    row=ifault/nstrike #Row number corresponding to this subfault
    column=ifault-(nstrike*row)
    if nstrike<4 or ndip<4:
        print "ERROR: The fault model is too small for Laplacian regualrization. You need a minimum of 4 rows and 4 columns in the model."
        return False,False
    if row==0 and column==0: #Top right corner
        stencil=array([ifault,ifault+1,ifault+nstrike])
        values=array([-4,1,1])
        if top.lower()=='free':
            values[2]+=1
        if right.lower=='free':
            values[1]+=1
        return stencil,values
    if row==0 and column==(nstrike-1): #Top left corner
        stencil=array([ifault,ifault-1,ifault+nstrike])
        values=array([-4,1,1])
        if top.lower()=='free':
            values[2]+=1
        if left.lower=='free':
            values[1]+=1
        return stencil,values
    if row==(ndip-1) and column==0: #Bottom right corner
        stencil=array([ifault,ifault+1,ifault-nstrike])
        values=([-4,1,1])
        if bottom.lower()=='free':
            values[2]+=1
        if right.lower=='free':
            values[1]+=1
        return stencil,values
    if row==(ndip-1) and column==(nstrike-1): #Bottom left corner
        stencil=array([ifault,ifault-1,ifault-nstrike])
        values=array([-4,1,1,])
        if bottom.lower()=='free':
            values[2]+=1
        if left.lower=='free':
            values[1]+=1
        return stencil,values
    if row==0: #Top edge, NOT the corner
        stencil=array([ifault,ifault+1,ifault-1,ifault+nstrike])
        values=array([-4,1,1,1])
        if top.lower()=='free':
            values[3]+=1
        return stencil,values
    if row==(ndip-1): #Bottom edge, NOT the corner
        stencil=array([ifault,ifault-1,ifault+1,ifault-nstrike])
        values=array([-4,1,1,1])
        if bottom.lower()=='free':
            values[3]+=1
        return stencil,values
    if column==0: #Right edge, NOT the corner
        stencil=array([ifault,ifault-nstrike,ifault+nstrike,ifault+1])
        values=array([-4,1,1,1])
        if right.lower()=='free':
            values[3]+=1
        return stencil,values
    if column==(nstrike-1): #left edge, NOT the corner
        stencil=array([ifault,ifault-nstrike,ifault+nstrike,ifault-1])
        values=array([-4,1,1,1])
        if left.lower()=='free':
            values[3]+=1
        return stencil,values
    else: #Somewhere in the middle
        stencil=array([ifault,ifault-1,ifault+1,ifault-nstrike,ifault+nstrike])
        values=array([-4,1,1,1,1])
        return stencil,values



def prep_synth(syn,st):
    '''
    Extend syntetic to start time of data and cut it to end time of data, make sure
    synthetic ALWAYS ends after data
    
    IN:
        syn: Synthetic stream object
        st: Data stream object
        
    OUT:
        syn: Trimmed synthetic
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
        
def gdims(datafiles,nfaults,decimate):
    '''
    Survey the data files to determine what dimension G will be and return a matrix of zeros 
    with the required dimensions
    '''
    
    from obspy import read
    from numpy import zeros
    from green import stdecimate
    
    npts=0
    if decimate==0:
        decimate=1
    for k in range(len(datafiles)):
        e=read(datafiles[k]+'.e')
        n=read(datafiles[k]+'.n')
        u=read(datafiles[k]+'.u')
        if e[0].stats.npts==n[0].stats.npts==u[0].stats.npts:
            if decimate>1:
                e=stdecimate(e,decimate)
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
    
    IN:
        outfile: Name of temporary station file
        sta: Single station you wish to include int he temp file
        lon: Station longitude
        lat: Station latitude
        gffiles: GF control file
        
    OUT:
        Nothing
    '''
    f=open(outfile,'a')
    for k in range(len(sta)):
        out=sta[k]+'\t'+repr(round(lon[k],6))+'\t'+repr(round(lat[k],6))+'\t'+gffiles[k]+'\n'
        f.write(out)
    f.close()

    
            
def epi2subfault(epicenter,source,vr,tr):
    '''
    Compute time delays from epicenter to subfault based on a give rupture speed Coordinates in 
    lat/lon,depth(km), vr in km/s, tr is delay to apply to rupture speed in secs.
    
    IN:
        epicenter: Epicentral coordinates
        source: Matrix of subfault coordinates
        vr: Rupture velocity
        tr: Timde delay to apply to rupture speed.
        
    OUT:
        tdelay: Time delays in seconds to all subfaults
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
    #Apply delay due to window
    tdelay=tdelay+tr
    return tdelay   
    
    

def get_stats(WG,sol,wd):
    '''
    Compute basic performance metrics of an inversion
    
    IN:
        WG: Dataweights times GFs WG=W*G
        sol: Soluction vector from inversion
        wd: Data weights times data vector, wd=W*d
    OUT:
        L2: ||Gm-d||
        Lm: ||Lm||
    '''
    
    from numpy.linalg import norm
    
    wds=WG.dot(sol)
    L2=norm(wds-wd)
    Lm=norm(sol)
    return L2,Lm
    
    
def get_VR(G,sol,d):
    '''
    Compute Variance reduction to the data
    
    IN:
        G: GF matrix
        sol: Solution  vector from inversion
        d: data vector
    OUT:
        VR: Variance reduction (%)
    '''
    
    ds=G.dot(sol)
    #Variance reduction
    res=((d-ds)**2)**0.5
    dnorm=(d**2)**0.5 #Yes i know this is dumb, shush
    VR=(1-(res.sum()/dnorm.sum()))*100
    return VR
    
    
def get_ABIC(G,GTG,sol,d,lambda_s,lambda_t,Ls,LsLs,Lt,LtLt):
    '''
    Compute Akaike's Bayesian information criterion, for details see Ide et al. (1996)
    in BSSA, specifically equation 33.
    
    IN:
        G: GFs matrix
        sol Solution vector from inversion
        d: Data vector
        lambda_s: Spatial regularization parameter
        lambda_t: Temporal regularization parameter
        Ls: Spatial regularization matrix
        Lt:Temporal reularization matrix
        Ls_rank: Rank of Ls (#eigenvalues>0)
        Lt_rank: Rank of Lt
    OUT:
        ABIC: Akaike's Bayesian information criterion
    
    '''
    
    from numpy import log
    from numpy.linalg  import norm,slogdet
    
    #Data points
    N=d.size
    #Model parameters
    M=sol.size
    #Off you go, compute it
    if lambda_t==0: #There is only one contraint (no temporal regularization)
        s=norm(d-G.dot(sol))**2+(lambda_s**2)*norm(Ls.dot(sol))**2
        a1=N*log(s)
        a2=M*log(lambda_s**2)
        sq,a3=slogdet(GTG+(lambda_s**2)*LsLs)
        #Add 'em up
        ABIC=a1-a2+a3
        return ABIC
    else: #There is a double regularization, use Fukahata et al. definition
        print '... computing 2d-ABIC'
        s=(norm(d-G.dot(sol))**2)+((lambda_s**2)*(norm(Ls.dot(sol))**2))+((lambda_t**2)*(norm(Lt.dot(sol))**2))
        a1=N*log(s)
        sq,a2=slogdet((lambda_s**2)*LsLs+(lambda_t**2)*LtLt)
        sq,a3=slogdet(GTG+(lambda_s**2)*LsLs+(lambda_t**2)*LtLt)
        #Add 'em up
        ABIC=a1-a2+a3
        return ABIC

    
    
def get_moment(home,project_name,fault_name,model_name,sol):
    '''
    Compute total moment from an inversion
    '''
    from numpy import log10,genfromtxt,loadtxt,arange,zeros
    from forward import get_mu
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Open structure file
    mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
    #Get slip quantities
    iss=2*arange(len(sol)/2)
    ids=2*arange(len(sol)/2)+1
    ss=sol[iss]
    ds=sol[ids]
    #Get total slip
    slip=(ss**2+ds**2)**0.5
    #Get rigidities and areas
    mu=zeros(len(ds))
    A=zeros(len(ds))
    i=0
    for krupt in range((len(sol)/len(f))/2):
        for k in range(len(f)):
            mu[i]=get_mu(mod,f[k,3])
            A[i]=f[k,8]*f[k,9]
            i+=1
    #Compute moments
    try:
        M0=mu*A*slip[:,0]
    except:
        M0=mu*A*slip
    #Total up and copute magnitude
    M0=M0.sum()
    Mw=(2./3)*(log10(M0)-9.1)
    return M0,Mw
    
    
def ds2rot(sol,beta):
    '''
    Rotate from a coordiante system where the basis is SS=0,DS=90 to one where the
    basis is SS=0+beta and DS=90+beta
    '''
    from numpy import array,deg2rad,cos,sin,arange,vstack,zeros
    
    #Split into strike-slip and dip-slip
    iss=arange(0,len(sol),2)
    ids=arange(1,len(sol),2)
    if len(iss)==1:
        ss=sol[0]
        ds=sol[1]
    else:
        ss=sol[iss]
        ds=sol[ids]
    #Rotate
    beta=deg2rad(beta)
    rot=array([[cos(beta),sin(beta)],[-sin(beta),cos(beta)]]).dot(vstack((ss,ds)))
    #Re-insert in output vector
    out=zeros(sol.shape)
    out[iss]=rot[0,:]
    out[ids]=rot[1,:]
    return out
    
def rot2ds(sol,beta):
    '''
    Reverses the operationd escribed in function ds2rot()
    '''
    from numpy import array,deg2rad,cos,sin,arange,vstack,zeros
    
    #Split into strike-slip and dip-slip
    iss=arange(0,len(sol),2)
    ids=arange(1,len(sol),2)
    if len(iss)==1:
        ssrot=sol[0]
        dsrot=sol[1]
    else:
        ssrot=sol[iss]
        dsrot=sol[ids]
    #Rotate
    beta=deg2rad(beta)
    ssds=array([[cos(beta),-sin(beta)],[sin(beta),cos(beta)]]).dot(vstack((ssrot.transpose(),dsrot.transpose())))
    #Re-insert in output vector
    out=zeros(sol.shape)
    out[iss,0]=ssds[0,:]
    out[ids,0]=ssds[1,:]
    return out
    
def resample_to_data(synth,data):
    '''
    Resample a synthetic to the time samples contained in the input data
    IN:
        synth: synthetic stream object
        data: data stream object
    OUT:
        st: resampled synthetic
    '''
    from datetime import timedelta
    from numpy import interp
    from forward import round_time
    
    #Get data sampling interval
    delta=data[0].stats.delta
    #Get synthetic start time
    t1synth=synth[0].stats.starttime
    #Get start time sampled at data interval and correction encessary
    t1=round_time(t1synth,delta)
    dt=t1-t1synth #This is the correctiont hat needs to be applied
    #Make current time vector
    tsynth=synth[0].times()
    tcorrect=tsynth+dt
    #Interpolate to new time vector
    synth_data=synth[0].data
    synth_correct=interp(tcorrect,tsynth,synth_data)
    #Place in output stream object
    synth[0].starttime=t1
    synth[0].data=synth_correct
    return synth