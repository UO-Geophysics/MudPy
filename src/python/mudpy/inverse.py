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
        num_windows,decimate,bandpass,tsunami=False,onset_file=None):
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
        decimate: Constant decimationf actor applied to GFs, set =0 for no decimation
        
    OUT:
        G: Fully assembled GF matrix
    '''
    
    from numpy import arange,genfromtxt,where,loadtxt,array,c_,concatenate,save,load,size,tile,expand_dims
    from os import remove
    from os.path import split
    
    G_name=home+project_name+'/GFs/matrices/'+G_name
    K_name=G_name+'.K'
    if G_from_file==True: #load from file
        if G_name[-3:]!='npy':
            K_name=K_name+'.npy'
            G_name=G_name+'.npy'
        print('Loading G from file '+G_name)
        G=load(G_name)
        #K=load(K_name)
    else: #assemble G one data type at a time
        print('Assembling G from synthetic computations...')
        
        #is there an onset times file?
        if onset_file is not None:
            tdelay_constant=genfromtxt(home+project_name+'/data/model_info/'+onset_file)
        
        #Read in GFlist and decide what to compute
        gf_file=home+project_name+'/data/station_info/'+GF_list
        mini_station=home+project_name+'/data/station_info/tempG.sta'
        stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='U')
        GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
        GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='U')
        #Check for single station sized arrays
        if GF.ndim==1: #One station
            GF=expand_dims(GF,0)
            GFfiles=expand_dims(GFfiles,0)
            stations=array([stations])
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
                #Gstatic=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay,decimate,lowpass)
                Ess=[] ; Eds=[] ; Nss=[] ; Nds=[] ; Zss=[] ; Zds=[]
                first_window=True
                Gstatic= makeG(home,project_name,fault_name,model_name,split(mini_station)[1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                remove(mini_station) #Cleanup  
        #Dispalcement waveform GFs
        kgf=3
        Gdisp=array([])
        if GF[:,kgf].sum()>0:
            #Load fault file
            source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
            trise=source[:,7]
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0 or len(array(i))>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,1])
                gftype='disp'
                #Decide on delays for each time window (50% overlap)
                delay_multiplier=tile(arange(0,num_windows),(len(trise),1)).T
                trupt=tile(trise/2,(num_windows,1))
                trupt=trupt*delay_multiplier
                for krup in range(num_windows):
                    print('Working on window '+str(krup+1))
                    if onset_file==None: #There is no onset times file, caluclate them
                        tdelay=epi2subfault(epicenter,source,rupture_speed,trupt[krup,:])
                    else: #use provided onset times and add time window delay
                        tdelay =  tdelay_constant+trupt[krup,:]
                    if krup==0: #First rupture speed
                        first_window=True
                        Ess=[] ; Eds=[] ; Nss=[] ; Nds=[] ; Zss=[] ; Zds=[]
                        Gdisp_temp,Ess,Eds,Nss,Nds,Zss,Zds = makeG(home,project_name,fault_name,model_name,split(mini_station)[1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                        Gdisp=Gdisp_temp
                    else:
                        first_window=False
                        Gdisp_temp,Ess,Eds,Nss,Nds,Zss,Zds = makeG(home,project_name,fault_name,model_name,split(mini_station)[1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                        Gdisp=c_[Gdisp,Gdisp_temp]
                remove(mini_station) #Cleanup 
        #Velocity waveforms
        kgf=4
        Gvel=array([])
        if GF[:,kgf].sum()>0:
            #Load fault file
            source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
            trise=source[:,7]
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0 or len(array(i))>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,2])
                gftype='vel'
                #Decide on delays for each time window (50% overlap)
                delay_multiplier=tile(arange(0,num_windows),(len(trise),1)).T
                trupt=tile(trise/2,(num_windows,1))
                trupt=trupt*delay_multiplier
                for krup in range(num_windows):
                    print('Working on window '+str(krup+1))
                    if onset_file==None: #There is no onset times file, caluclate them
                        tdelay=epi2subfault(epicenter,source,rupture_speed,trupt[krup,:])
                    else: #use provided onset times and add time window delay
                        tdelay =  tdelay_constant+trupt[krup,:]
                    if krup==0: #First rupture speed
                        first_window=True
                        Ess=[] ; Eds=[] ; Nss=[] ; Nds=[] ; Zss=[] ; Zds=[]
                        Gvel_temp,Ess,Eds,Nss,Nds,Zss,Zds = makeG(home,project_name,fault_name,model_name,mini_station.split('/')[-1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                        Gvel=Gvel_temp
                    else:
                        first_window=False
                        Gvel_temp,Ess,Eds,Nss,Nds,Zss,Zds = makeG(home,project_name,fault_name,model_name,mini_station.split('/')[-1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                        Gvel=c_[Gvel,Gvel_temp]
                remove(mini_station) #Cleanup 
        #Tsunami waveforms
        kgf=5
        Gtsun=array([])
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
            if len(i)>0 or len(array(i))>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,3])
                gftype='tsun'
                #Decide on delays for each time window (50% overlap)
                trupt=arange(0,num_windows)*trise/2
                for krup in range(num_windows):
                    print('Working on window '+str(krup+1))
                    if onset_file==None: #There is no onset times file, caluclate them
                        tdelay=epi2subfault(epicenter,source,rupture_speed,trupt[krup])
                    else: #use provided onset times and add time window delay
                        tdelay =  tdelay_constant+trupt[krup]
                    if krup==0: #First rupture speed
                        first_window=True
                        Ess=[] ; Eds=[] ; Nss=[] ; Nds=[] ; Zss=[] ; Zds=[]
                        Gtsun_temp,SS,DS = makeG(home,project_name,fault_name,model_name,mini_station.split('/')[-1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                        Gtsun=Gtsun_temp
                    else:
                        first_window=False
                        Gtsun_temp,SS,DS = makeG(home,project_name,fault_name,model_name,split(mini_station)[1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,SS,DS,Nss,Nds,Zss,Zds)
                        Gtsun=c_[Gtsun,Gtsun_temp]
                remove(mini_station) #Cleanup 
        #InSAR LOS offsets
        kgf=6
        Ginsar=array([])
        if GF[:,kgf].sum()>0:
            try:
                remove(mini_station) #Cleanup  
            except:
                pass
            #Make mini station file 
            i=where(GF[:,kgf]==1)[0]
            if len(i)>0: #There's actually something to do
                mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,0])
                gftype='insar'
                tdelay=0
                #Gstatic=makeG(home,project_name,fault_name,model_name,split(mini_station)[1],gftype,tdelay,decimate,lowpass)
                Ess=[] ; Eds=[] ; Nss=[] ; Nds=[] ; Zss=[] ; Zds=[]
                first_window=True
                Ginsar= makeG(home,project_name,fault_name,model_name,split(mini_station)[1],
                                                                gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds)
                remove(mini_station) #Cleanup     
            
        #Done, now concatenate them all ccompute transpose product and save
        if num_windows==1: #Only one window
            G=concatenate([g for g in [Gstatic,Gdisp,Gvel,Gtsun,Ginsar] if g.size > 0])
        elif num_windows>1: #Multiple windows, this means we need tot ile the statics if they exist
            if size(Gstatic)!=0: #Static is not empty, need to tile it
                Gstatic_nwin=Gstatic.copy()
                for nwin in range(num_windows-1):
                    Gstatic_nwin=c_[Gstatic_nwin,Gstatic]
                Gstatic=Gstatic_nwin.copy()
                Gstatic_nwin=None #Release memory
            if size(Ginsar)!=0: #Static is not empty, need to tile it
                Ginsar_nwin=Ginsar.copy()
                for nwin in range(num_windows-1):
                    Ginsar_nwin=c_[Ginsar_nwin,Ginsar]
                Ginsar=Ginsar_nwin.copy()
                Ginsar_nwin=None #Release memory
            print(Gstatic.shape)
            print(Gvel.shape)
            print(Ginsar.shape)
            G=concatenate([g for g in [Gstatic,Gdisp,Gvel,Gtsun,Ginsar] if g.size > 0])
        print('Saving GF matrix to '+G_name+' this might take just a second...')
        save(G_name,G)
        #save(K_name,K)
    return G
    
    
def makeG(home,project_name,fault_name,model_name,station_file,gftype,tsunami,tdelay,decimate,bandpass,first_window,Ess,Eds,Nss,Nds,Zss,Zds):
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
    from numpy import genfromtxt,loadtxt,zeros,array,inf,size,where
    from obspy import read,Stream,Trace
    from mudpy.forward import tshift
    from mudpy.forward import lowpass as lfilt
    from mudpy.forward import highpass as hfilt
    from mudpy.green import stdecimate
    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    datafiles=genfromtxt(station_file,dtype="U",usecols=3)
    try:
        Nsta=len(staname)
    except:
        Nsta=1
        staname=array([staname])
        datafiles=array([datafiles])
    insert_position=0
    #Initalize G for faster assignments
    if gftype.lower()=='static': #Initialize output matrix
        G=zeros((Nsta*3,Nfaults*2))
    elif gftype.lower()=='insar': #Initialize output matrix
        G=zeros((Nsta,Nfaults*2))
    else:
        pass #For disp or vel waveforms G is initalized below
    if gftype.lower()=='static': #Make matrix of static GFs
        for ksta in range(Nsta):
            print('Assembling static GFs for station '+staname[ksta])
            #Initalize output variable
            Gtemp=zeros([3,Nfaults*2])
            #Where's the data
            statics_path=home+project_name+'/GFs/static/'
            #Loop over subfaults
            for kfault in range(Nfaults):
                
                nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                strdepth='%.4f' % source[kfault,3]
                syn_path=statics_path+model_name+'_'+strdepth+'.'+nsub+'/'
                
                                
                if kfault%10==0:
                    print('... working on subfault '+str(kfault)+' of '+str(Nfaults))
                
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                
                coseis_ss=loadtxt(syn_path+nfault+'.SS.static.neu',usecols=[1,2,3])
                coseis_ds=loadtxt(syn_path+nfault+'.DS.static.neu',usecols=[1,2,3])
                coseis_sta=loadtxt(syn_path+nfault+'.SS.static.neu',usecols=0,dtype='U')
                                
                #Find row corresponding to station
                ista=where(coseis_sta==staname[ksta])[0]
                
                #get values
                nss=coseis_ss[ista,0]
                ess=coseis_ss[ista,1]
                zss=coseis_ss[ista,2]
                
                nds=coseis_ds[ista,0]
                eds=coseis_ds[ista,1]
                zds=coseis_ds[ista,2]

                
                #Place into G matrix
                Gtemp[0,2*kfault]=nss   ; Gtemp[0,2*kfault+1]=nds    #North
                Gtemp[1,2*kfault]=ess ; Gtemp[1,2*kfault+1]=eds  #East
                Gtemp[2,2*kfault]=zss ; Gtemp[2,2*kfault+1]=zds  #Up
                #Append to G
            #Append to output matrix
            G[ksta*3:ksta*3+3,:]=Gtemp   
        return G   
    if gftype.lower()=='insar': #Make matrix of insar LOS GFs
        for ksta in range(Nsta):
            print('Assembling static GFs for station '+staname[ksta])
            #Initalize output variable
            Gtemp=zeros([1,Nfaults*2])
            #Where's the data
            statics_path=home+project_name+'/GFs/static/'
            #Data path, need this to find LOS vector
            los_path=home+project_name+'/data/statics/'
            #Read los vector for this subfault
            los=genfromtxt(los_path+staname[ksta]+'.los')
            los=los[1:]
            #Loop over subfaults
            for kfault in range(Nfaults):
                
                nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                strdepth='%.4f' % source[kfault,3]
                syn_path=statics_path+model_name+'_'+strdepth+'.'+nsub+'/'
                
                
                if kfault%10==0:
                    print('... working on subfault '+str(kfault)+' of '+str(Nfaults))
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                
                coseis_ss=loadtxt(syn_path+nfault+'.SS.static.neu',usecols=[1,2,3])
                coseis_ds=loadtxt(syn_path+nfault+'.DS.static.neu',usecols=[1,2,3])
                coseis_sta=loadtxt(syn_path+nfault+'.SS.static.neu',usecols=0,dtype='U')
                
                #Find row corresponding to station
                ista=where(coseis_sta==staname[ksta])[0]
                
                #get values
                nss=coseis_ss[ista,0]
                ess=coseis_ss[ista,1]
                zss=coseis_ss[ista,2]
                
                nds=coseis_ds[ista,0]
                eds=coseis_ds[ista,1]
                zds=coseis_ds[ista,2]
                
                # Dot product of GFs and los vector
                los_ss=los.dot(array([nss,ess,zss]))
                los_ds=los.dot(array([nds,eds,zds]))
                
                #Place into G matrix
                Gtemp[0,2*kfault]=los_ss   ; Gtemp[0,2*kfault+1]=los_ds  
                #Append to G
            #Append to output matrix
            G[ksta*1,:]=Gtemp   
        return G  
    if gftype.lower()=='disp' or gftype.lower()=='vel':  #Full waveforms
        if gftype.lower()=='disp':
            vord='disp'
            if bandpass == None:
                BP=None
            else:
                BP=bandpass[0]
        else:
            vord='vel'
            if bandpass == None:
                BP=None
            else:
                BP=bandpass[1]
            tsunami=False

    #Tsunami and requires bandpassing?
    if gftype.lower()=='disp' or gftype.lower()=='vel':  #Full waveforms
        if gftype.lower()=='disp':
            vord='disp'
            if bandpass == None:
                BP=None
            else:
                BP=bandpass[0]

        if first_window==True: #Read in GFs from file
            ktrace=0
            for ksta in range(Nsta):
                print('Reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
                for kfault in range(Nfaults):
                    #Get subfault GF directory
                    nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                    nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                    strdepth='%.4f' % source[kfault,3]
                    if tsunami==True:
                        syn_path=home+project_name+'/GFs/tsunami/'+model_name+'_'+strdepth+'.'+nsub+'/'
                    else:
                        syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                    #Get synthetics
                    if kfault==0 and ksta==0: #It's the first one, initalize stream object
                        Ess=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
                        Nss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
                        Zss=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
                        Eds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
                        Nds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
                        Zds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
                    else: #Just add to stream object
                        Ess+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
                        Nss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
                        Zss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
                        Eds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
                        Nds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
                        Zds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
                    #Perform operations that need to only happen once (filtering and decimation)
                    if any(BP!=None):# or ksta==1: #Apply filter to GFs
                        if kfault==0:
                            print('Bandpassing on frequency band:')
                            print('... '+str(BP))
                            print('... station '+staname[ksta])
                        fsample=1./Ess[ktrace].stats.delta
                        if len(BP)>1 and BP[1]==inf: #A high pass filter has been requested
                            if kfault==0:
                                print('.... highpass')
                            Ess[ktrace].data=hfilt(Ess[ktrace].data,BP[0],fsample,2)
                            Nss[ktrace].data=hfilt(Nss[ktrace].data,BP[0],fsample,2)
                            Zss[ktrace].data=hfilt(Zss[ktrace].data,BP[0],fsample,2)
                            Eds[ktrace].data=hfilt(Eds[ktrace].data,BP[0],fsample,2)
                            Nds[ktrace].data=hfilt(Nds[ktrace].data,BP[0],fsample,2)
                            Zds[ktrace].data=hfilt(Zds[ktrace].data,BP[0],fsample,2)                        
                        else: #A low pass or bandpass filter has been requested
                            if kfault==0:
                                print('.... lowpass')
                            #Ess[ktrace].data=lfilt(Ess[ktrace].data,BP[0],fsample,2)
                            #Nss[ktrace].data=lfilt(Nss[ktrace].data,BP[0],fsample,2)
                            #Zss[ktrace].data=lfilt(Zss[ktrace].data,BP[0],fsample,2)
                            #Eds[ktrace].data=lfilt(Eds[ktrace].data,BP[0],fsample,2)
                            #Nds[ktrace].data=lfilt(Nds[ktrace].data,BP[0],fsample,2)
                            #Zds[ktrace].data=lfilt(Zds[ktrace].data,BP[0],fsample,2)
                            Ess[ktrace].data=lfilt(Ess[ktrace].data,BP,fsample,2)
                            Nss[ktrace].data=lfilt(Nss[ktrace].data,BP,fsample,2)
                            Zss[ktrace].data=lfilt(Zss[ktrace].data,BP,fsample,2)
                            Eds[ktrace].data=lfilt(Eds[ktrace].data,BP,fsample,2)
                            Nds[ktrace].data=lfilt(Nds[ktrace].data,BP,fsample,2)
                            Zds[ktrace].data=lfilt(Zds[ktrace].data,BP,fsample,2)
                        #bandpass=None
                    
                    ### HAAAAAACCCCKKKK!!!!
                    #if vord=='vel': #Apply low pass filter to data (** DIRTY HACK!**)
                    #    if kfault==0:
                    #        print(vord
                    #        print('Bandpassing hack on velocity...'
                    #    bandpass=[1./20,1./2]
                    #    fsample=1./Ess[0].stats.delta
                    #    Ess[ktrace].data=lfilt(Ess[ktrace].data,bandpass,fsample,2)
                    #    Nss[ktrace].data=lfilt(Nss[ktrace].data,bandpass,fsample,2)
                    #    Zss[ktrace].data=lfilt(Zss[ktrace].data,bandpass,fsample,2)
                    #    Eds[ktrace].data=lfilt(Eds[ktrace].data,bandpass,fsample,2)
                    #    Nds[ktrace].data=lfilt(Nds[ktrace].data,bandpass,fsample,2)
                    #    Zds[ktrace].data=lfilt(Zds[ktrace].data,bandpass,fsample,2)
                    #    #bandpass=None
                    ### END AWFUL TERRIBLE REALLY VERY BAD HACK
                        
                    if decimate!=None: 
                        Ess[ktrace]=stdecimate(Ess[ktrace],decimate)
                        Nss[ktrace]=stdecimate(Nss[ktrace],decimate)
                        Zss[ktrace]=stdecimate(Zss[ktrace],decimate)    
                        Eds[ktrace]=stdecimate(Eds[ktrace],decimate)
                        Nds[ktrace]=stdecimate(Nds[ktrace],decimate)
                        Zds[ktrace]=stdecimate(Zds[ktrace],decimate)
                        
                    
                    ktrace+=1            
        #Read time series
        for ksta in range(Nsta):
            edata=read(datafiles[ksta]+'.e')
            ndata=read(datafiles[ksta]+'.n')
            udata=read(datafiles[ksta]+'.u')
            if decimate!=None:
                edata[0]=stdecimate(edata[0],decimate)
                ndata[0]=stdecimate(ndata[0],decimate)
                udata[0]=stdecimate(udata[0],decimate)
            if ksta==0:
                Edata=edata.copy()
                Ndata=ndata.copy()
                Udata=udata.copy()
            else:
                Edata+=edata
                Ndata+=ndata
                Udata+=udata            
        #Finished reading, filtering, etc now time shift by rupture time and resmaple to data
        ktrace=0
        print("Aligning GFs and resampling to data times...")
        for ksta in range(Nsta):
            #Loop over subfaults
            print('...Working on station #'+str(ksta+1)+' of '+str(Nsta))
            for kfault in range(Nfaults):
                #Assign current GFs
                ess=Stream(Trace())
                nss=Stream(Trace())
                zss=Stream(Trace())
                eds=Stream(Trace())
                nds=Stream(Trace())
                zds=Stream(Trace())
                ess[0]=Ess[ktrace].copy()
                nss[0]=Nss[ktrace].copy()
                zss[0]=Zss[ktrace].copy()
                eds[0]=Eds[ktrace].copy()
                nds[0]=Nds[ktrace].copy()
                zds[0]=Zds[ktrace].copy()
                #Time shift them according to subfault rupture time, zero pad, round to dt interval,decimate
                #and extend to maximum time
                ess=tshift(ess,tdelay[kfault])
                nss=tshift(nss,tdelay[kfault])
                zss=tshift(zss,tdelay[kfault])
                eds=tshift(eds,tdelay[kfault])
                nds=tshift(nds,tdelay[kfault])
                zds=tshift(zds,tdelay[kfault])
                #Now time align stuff                                
                ess=resample_to_data(ess[0],Edata[ksta])
                ess=prep_synth(ess,Edata[ksta])
                nss=resample_to_data(nss[0],Ndata[ksta])
                nss=prep_synth(nss,Ndata[ksta])
                zss=resample_to_data(zss[0],Udata[ksta])
                zss=prep_synth(zss,Udata[ksta])
                eds=resample_to_data(eds[0],Edata[ksta])
                eds=prep_synth(eds,Edata[ksta])
                nds=resample_to_data(nds[0],Ndata[ksta])
                nds=prep_synth(nds,Ndata[ksta])
                zds=resample_to_data(zds[0],Udata[ksta])
                zds=prep_synth(zds,Udata[ksta])
                #Insert into Gtemp then append to G
                if kfault==0 and ksta==0: #It's the first subfault and station, initalize G
                    G=gdims(datafiles,Nfaults,decimate) #Survey all stations to decide size of G
                if kfault==0: #Initalize Gtemp (different size for each station)
                    #How many points left in the tiem series
                    npts=Edata[ksta].stats.npts
                    print("... ... "+str(npts)+" data points left over after decimation")
                    Gtemp=zeros([3*npts,Nfaults*2])      
                #Insert synthetics into Gtemp
                Gtemp[0:npts,2*kfault]=nss.data
                Gtemp[0:npts,2*kfault+1]=nds.data
                Gtemp[npts:2*npts,2*kfault]=ess.data
                Gtemp[npts:2*npts,2*kfault+1]=eds.data
                Gtemp[2*npts:3*npts,2*kfault]=zss.data
                Gtemp[2*npts:3*npts,2*kfault+1]=zds.data
                ktrace+=1
            #After looping through all subfaults Insert Gtemp into G
            G[insert_position:insert_position+3*npts,:]=Gtemp
            insert_position+=3*npts #Update for next station
        return G,Ess,Eds,Nss,Nds,Zss,Zds
    
    if gftype.lower()=='tsun':
        if first_window==True: #Read in GFs from file
            ktrace=0
            for ksta in range(Nsta):
                print('Reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
                for kfault in range(Nfaults):
                    #Get subfault GF directory
                    nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                    nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                    strdepth='%.4f' % source[kfault,3]
                    syn_path=home+project_name+'/GFs/tsunami/'+model_name+'_'+strdepth+'.'+nsub+'/'
                    #Get synthetics
                    if kfault==0 and ksta==0: #It's the first one, initalize stream object
                        SS=read(syn_path+staname[ksta]+'.ss.tsun')
                        DS=read(syn_path+staname[ksta]+'.ds.tsun')
                    else: #Just add to stream object
                        SS+=read(syn_path+staname[ksta]+'.ss.tsun')
                        DS+=read(syn_path+staname[ksta]+'.ds.tsun')
                    
                    #Band pass filter data
                    if bandpass == None:
                        BP=None
                    else:
                        BP=bandpass[2]
                    
                    if size(BP)>1:
                        if kfault==0:
                            print('Bandpassing on frequency band:')
                            print('... '+str(BP))
                            print('... station '+staname[ksta])
                        fsample=1./SS[ktrace].stats.delta
                        if BP[1]==inf: #A high pass filter has been requested
                            SS[ktrace].data=hfilt(SS[ktrace].data,BP[0],fsample,2)
                            DS[ktrace].data=hfilt(DS[ktrace].data,BP[0],fsample,2)                    
                        else: #A low pass or bandpass filter has been requested
                            SS[ktrace].data=lfilt(SS[ktrace].data,BP,fsample,2)
                            DS[ktrace].data=lfilt(DS[ktrace].data,BP,fsample,2)
                                                
                    ktrace+=1            
        else:
            SS=Ess.copy()
            DS=Eds.copy()
        #Read time series
        for ksta in range(Nsta):
            if ksta==0:
                Data=read(datafiles[ksta])
            else:
                Data+=read(datafiles[ksta])         
        #Finished reading, filtering, etc now time shift by rupture time and resmaple to data
        ktrace=0
        print("Aligning GFs and resampling to data times...")
        for ksta in range(Nsta):
            #Loop over subfaults
            print('...Working on station #'+str(ksta+1)+' of '+str(Nsta))
            for kfault in range(Nfaults):
                #Assign current GFs
                ss=Stream(Trace())
                ds=Stream(Trace())
                ss[0]=SS[ktrace].copy()
                ds[0]=DS[ktrace].copy()
                #Time shift them according to subfault rupture time, zero pad, round to dt interval,decimate
                #and extend to maximum time
                ss=tshift(ss,tdelay[kfault])
                ds=tshift(ds,tdelay[kfault])
                #Now time align stuff                                
                #ss=resample_synth_tsun(ss[0],Data[ksta])
                ss=prep_synth(ss[0],Data[ksta])
                #ds=resample_synth_tsun(ds[0],Data[ksta])
                ds=prep_synth(ds[0],Data[ksta])
                #Insert into Gtemp then append to G
                if kfault==0 and ksta==0: #It's the first subfault and station, initalize G
                    G=gdims_tsun(datafiles,Nfaults,decimate) #Survey all stations to decide size of G
                if kfault==0: #Initalize Gtemp (different size for each station)
                    #How many points left in the tiem series
                    npts=Data[ksta].stats.npts
                    print("... ... "+str(npts)+" data points left over")
                    Gtemp=zeros([npts,Nfaults*2])      
                #Insert synthetics into Gtempview
                Gtemp[0:npts,2*kfault]=ss.data
                Gtemp[0:npts,2*kfault+1]=ds.data
                ktrace+=1
            #After looping through all subfaults Insert Gtemp into G
            G[insert_position:insert_position+npts,:]=Gtemp
            insert_position+=npts #Update for next station
        return G,SS,DS             
    if gftype.lower()=='strain':
        pass                                
      
      
def getdata(home,project_name,GF_list,decimate,bandpass,quiet=False):
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
    from mudpy.forward import round_time
    from mudpy.green import stdecimate
    from mudpy.forward import lowpass as lfilt
    

    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='U')  
    stations=genfromtxt(gf_file,usecols=0,dtype='U')  
    
    #Parse out filtering pass bands
    displacement_bandpass=bandpass[0]
    velocity_bandpass=bandpass[1]
    tsunami_bandpass=bandpass[2]
    
    #Read one column at a time
    kgf=0 #Static field
    dstatic=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        if quiet==False:  
            print('Assembling static offsets from '+stations[i[ksta]]+' into data vector.')
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
        if quiet==False:  
            print('Assembling displacement waveforms from '+str(stations[i[ksta]])+' into data vector.')
        #print(str(GFfiles[i[ksta],kgf]))
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        if displacement_bandpass is not None: #Apply low pass filter to data
            fsample=1./n[0].stats.delta
            n[0].data=lfilt(n[0].data,displacement_bandpass,fsample,2)
            e[0].data=lfilt(e[0].data,displacement_bandpass,fsample,2)
            u[0].data=lfilt(u[0].data,displacement_bandpass,fsample,2)
        #Decimate
        if decimate!=None:
            n[0]=stdecimate(n[0],decimate)
            e[0]=stdecimate(e[0],decimate)
            u[0]=stdecimate(u[0],decimate)
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
        if quiet==False:  
            print('Assembling acceleration waveforms from '+stations[i[ksta]]+' into data vector.')
        n=read(GFfiles[i[ksta],kgf]+'.n')
        e=read(GFfiles[i[ksta],kgf]+'.e')
        u=read(GFfiles[i[ksta],kgf]+'.u')
        if velocity_bandpass is not None: #Apply low pass filter to data
            fsample=1./n[0].stats.delta
            n[0].data=lfilt(n[0].data,velocity_bandpass,fsample,2)
            e[0].data=lfilt(e[0].data,velocity_bandpass,fsample,2)
            u[0].data=lfilt(u[0].data,velocity_bandpass,fsample,2)
        #Decimate
        if decimate!=None:
            n[0]=stdecimate(n[0],decimate)
            e[0]=stdecimate(e[0],decimate)
            u[0]=stdecimate(u[0],decimate)
        #Make sure they are rounded to a dt interval
        dt=e[0].stats.delta
        e[0].stats.starttime=round_time(e[0].stats.starttime,dt)
        n[0].stats.starttime=round_time(n[0].stats.starttime,dt)
        u[0].stats.starttime=round_time(u[0].stats.starttime,dt)
        dvel=append(dvel,r_[n[0].data,e[0].data,u[0].data])
    
    #Tsunami
    kgf=3
    dtsun=array([])
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        if quiet==False:  
            print('Assembling tsunami waveforms from '+stations[i[ksta]]+' into data vector.')
        tsun=read(GFfiles[i[ksta],kgf])
        if tsunami_bandpass is not None: #Apply low pass filter to data
            fsample=1./tsun[0].stats.delta
            tsun[0].data=lfilt(tsun[0].data,tsunami_bandpass,fsample,2)
            
        dtsun=append(dtsun,tsun)
    
    #InSAR LOS
    kgf=4
    dlos=array([])  
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        if quiet==False:  
            print('Assembling InSAR LOS offsets from '+stations[i[ksta]]+' into data vector.')
        dtemp=genfromtxt(GFfiles[i[ksta],kgf])
        los=dtemp[0]
        dlos=append(dlos,los)          
    #Done, concatenate all, convert to column vector and exit
    d=concatenate([dx for dx in [dstatic,ddisp,dvel,dtsun,dlos] if dx.size > 0])
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
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    N=len(source) #No. of subfaults
    nstrike=nfaults[0]
    ndip=nfaults[1]
    #Initalize
    L=zeros((2*N,2*N))
    #Which L am I building?
    print('Making discrete Laplace operator regularization matrix...')
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
        print('WARNING temporal regularization is unecessary when only employing 1 time window. Returning zeros.')
        Lout=zeros((2*N,2*N))
        return Lout
    #Initalize
    L=eye(2*N*num_windows)
    #Ltest=zeros((N*num_windows,2*N*num_windows))
    print('Making first derivative temporal regularization matrix...')
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
    from mudpy.green import stdecimate

    print('Computing data weights...')
    #Read gf file and decide what needs tog et loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='U')
    weights=genfromtxt(gf_file,usecols=range(13,28),dtype='f')
    #Initalize
    w=zeros(len(d))
    kinsert=0
    #Deal with decimation
    #if decimate==None:
    #    decimate=1
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
        if decimate!=None:
            st[0]=stdecimate(st[0],decimate)
        nsamples=st[0].stats.npts
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
        if decimate!=None:
            st[0]=stdecimate(st[0],decimate)
        nsamples=st[0].stats.npts
        wn=(1/weights[i[ksta],6])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=wn
        kinsert=kinsert+nsamples
        we=(1/weights[i[ksta],7])*ones(nsamples)
        w[kinsert:kinsert+nsamples]=we
        kinsert=kinsert+nsamples
        wu=(1/weights[i[ksta],8])*ones(nsamples)
        try:
            w[kinsert:kinsert+nsamples]=wu
        except:
            print('WARNING: mismatch in data weights length')
            w[kinsert:kinsert+nsamples]=wu[:-1]  #I don't know why this happens, wu is the right length but w is too short by one but only for the last station
        kinsert=kinsert+nsamples
    #Tsunami
    kgf=3
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        #Read waveform to determine length of insert
        st=read(GFfiles[i[ksta],kgf])
        nsamples=st[0].stats.npts
        wtsun=(1/weights[i[ksta],9])*ones(nsamples)
        
#        wtsun2=read('/Users/dmelgarm/Oaxaca2020/tsunami/dart/43413.weights.sac')
#        wtsun2=1/wtsun2[0].data 
#        print('WARNING: Tsunami weights have been hard coded')
#        w[kinsert:kinsert+nsamples]=wtsun*wtsun2
        
        w[kinsert:kinsert+nsamples]=wtsun
        kinsert=kinsert+nsamples
    #InSAR
    kgf=4
    i=where(GF[:,kgf]==1)[0]
    for ksta in range(len(i)):
        w[kinsert]=1/weights[i[ksta],10] #LOS
        kinsert+=1
    #Make W and exit
    return w

    
    
#=================        Write inversion results      =========================
    
def write_model(home,project_name,run_name,fault_name,model_name,rupture_speed,num_windows,epicenter,sol,num,onset_file=None):
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
    from mudpy.forward import get_mu
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    trise=f[0,7]
    #Open structure file
    mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
    
    # are there forced onset times?
    if onset_file is not None:
        tdelay_constant = genfromtxt(home+project_name+'/data/model_info/'+onset_file)
    
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
        if onset_file==None: #No forced onset times
            trup[krup*(len(ds)//num_windows):(krup+1)*(len(ds)//num_windows)]=epi2subfault(epicenter,f,rupture_speed,trupt[krup])
        else: # use specified onset times plus window time delay
            trup[krup*(len(ds)//num_windows):(krup+1)*(len(ds)//num_windows)]=tdelay_constant+trupt[krup]
    #Prepare for output
    out1=f[:,0:8]
    out2=f[:,8:10]
    for k in range(num_windows-1):
        out1=r_[out1,f[:,0:8]]
        out2=r_[out2,f[:,8:10]]
    out=c_[out1,ss,ds,out2,trup,mu]
    outdir=home+project_name+'/output/inverse_models/models/'+run_name+'.'+str(num).rjust(4,'0')+'.inv'
    #CHANGE this to rupture definition as #No  x            y        z(km)      str     dip      rake       rise    dura     slip    ss_len  ds_len rupt_time
    fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e\t%10.1f\t%10.1f\t%8.4f\t%.4e'
    print('... writing model results to file '+outdir)
    savetxt(outdir,out,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
        
def write_synthetics_GOD(home,project_name,run_name,GF_list,ds,num,decimate):
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
    from numpy import array,savetxt,genfromtxt,squeeze,r_
    from mudpy.green import stdecimate
    
    print('... computing and saving synthetics...')
    num=str(num).rjust(4,'0')
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],skip_header=1,dtype='U')
    #Separate into its constituent parts (statics,displacaments, velocities, etc...)
    kinsert=0
#    kgf=0
    for i in range(len(stations)):
        sta=stations[i]
        if GF[i][0]==1: # STATIC DISPLACEMENT DATA
            neu=array([ds[kinsert],ds[kinsert+1],ds[kinsert+2]])
            #print "... ... " + sta + " static synthetics saved: " + str(kinsert) + " - " + str(kinsert+3) 
            kinsert+=3
            savetxt(home+project_name+'/output/inverse_models/statics/'+run_name+'.'+num+'.'+sta+'.static.neu',neu)
        if GF[i][1]==1: # DISPLACEMENT WAVEFORM DATA
            n=read(GFfiles[i,1]+'.n')
            e=read(GFfiles[i,1]+'.e')
            u=read(GFfiles[i,1]+'.u')
            if decimate != None:
                n[0]=stdecimate(n[0],decimate)
                e[0]=stdecimate(e[0],decimate)
                u[0]=stdecimate(u[0],decimate)
            npts=n[0].stats.npts
            n[0].data=squeeze(ds[kinsert:kinsert+npts])
            e[0].data=squeeze(ds[kinsert+npts:kinsert+2*npts])
            u[0].data=squeeze(ds[kinsert+2*npts:kinsert+3*npts])
            #print "... ... " + sta + " displacement synthetics saved: " + str(kinsert) + " - " + str(kinsert+3*npts)
            kinsert+=3*npts
            n.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.n.sac',format='SAC')
            e.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.e.sac',format='SAC')
            u.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.disp.u.sac',format='SAC')
        if GF[i][2]==1: #VELOCITY WAVEFORM DATA
            n=read(GFfiles[i,2]+'.n')
            e=read(GFfiles[i,2]+'.e')
            u=read(GFfiles[i,2]+'.u')
            if decimate != None:
                n[0]=stdecimate(n[0],decimate)
                e[0]=stdecimate(e[0],decimate)
                u[0]=stdecimate(u[0],decimate)
            npts=n[0].stats.npts
            n[0].data=squeeze(ds[kinsert:kinsert+npts])
            e[0].data=squeeze(ds[kinsert+npts:kinsert+2*npts])
            u[0].data=squeeze(ds[kinsert+2*npts:kinsert+3*npts])
            #print "... ... " + sta + " velocity synthetics saved: " + str(kinsert) + " - " + str(kinsert+3*npts)
            kinsert+=3*npts
            n.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.n.sac',format='SAC')
            e.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.e.sac',format='SAC')
            u.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.u.sac',format='SAC')
        if GF[i][3]==1: #TSUNAMI DATA
            tsun=read(GFfiles[i,3])
            npts=tsun[0].stats.npts
            synth=tsun.copy()
            synth[0].data=squeeze(ds[kinsert:kinsert+npts])
            #print "... ... " + sta + " tsunami synthetics saved: " + str(kinsert) + " - " + str(kinsert+npts)
            kinsert+=npts
            synth.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.tsun',format='SAC')
        if GF[i][4]==1: #INSAR DATA
            los_vector=genfromtxt(GFfiles[i,4])
            los_vector=los_vector[1:]
            los=array(r_[ds[kinsert],los_vector])
            #print "... ... " + sta + " InSAR synthetics saved: " + str(kinsert) + " - " + str(kinsert+1)
            kinsert+=1
            savetxt(home+project_name+'/output/inverse_models/statics/'+run_name+'.'+num+'.'+sta+'.los',los)

    
def write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,num,decimate):
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
    from numpy import array,savetxt,where,genfromtxt,squeeze,r_
    from mudpy.green import stdecimate
    
    print('... computing and saving synthetics...')
    num=str(num).rjust(4,'0')
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],skip_header=1,dtype='U')
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
            if decimate != None:
                n[0]=stdecimate(n[0],decimate)
                e[0]=stdecimate(e[0],decimate)
                u[0]=stdecimate(u[0],decimate)
            npts=n[0].stats.npts
            n[0].data=squeeze(ds[kinsert:kinsert+npts])
            e[0].data=squeeze(ds[kinsert+npts:kinsert+2*npts])
            u[0].data=squeeze(ds[kinsert+2*npts:kinsert+3*npts])
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
            if decimate !=None:
                n[0]=stdecimate(n[0],decimate)
                e[0]=stdecimate(e[0],decimate)
                u[0]=stdecimate(u[0],decimate)
            npts=n[0].stats.npts
            n[0].data=squeeze(ds[kinsert:kinsert+npts])
            e[0].data=squeeze(ds[kinsert+npts:kinsert+2*npts])
            u[0].data=squeeze(ds[kinsert+2*npts:kinsert+3*npts])
            kinsert+=3*npts
            n.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.n.sac',format='SAC')
            e.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.e.sac',format='SAC')
            u.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.vel.u.sac',format='SAC')
    #Tsunami
    kgf=3
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            tsun=read(GFfiles[i[ksta],kgf])
            npts=tsun[0].stats.npts
            synth=tsun.copy()
            synth[0].data=squeeze(ds[kinsert:kinsert+npts])
            kinsert+=npts
            synth.write(home+project_name+'/output/inverse_models/waveforms/'+run_name+'.'+num+'.'+sta+'.tsun',format='SAC')
    # InSAR
    kgf=4
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            los_vector=genfromtxt(home+project_name+'/data/statics/'+stations[i[ksta]]+'.los')
            los_vector=los_vector[1:]
            los=array(r_[ds[kinsert],los_vector])
            kinsert+=1
            savetxt(home+project_name+'/output/inverse_models/statics/'+run_name+'.'+num+'.'+sta+'.los',los)
            
        
def write_log(home,project_name,run_name,k,rupture_speed,num_windows,lambda_spatial,lambda_temporal,
        beta,L2,Lm,VR,ABIC,Mo,Mw,velmod,fault,g_name,gflist,solver,L2data):
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
    
    
    num=str(k).rjust(4,'0')
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
    f.write('VR static(%) = '+repr(VR[0])+'\n')
    f.write('VR displacement(%) = '+repr(VR[1])+'\n')
    f.write('VR velocity(%) = '+repr(VR[2])+'\n')
    f.write('VR tsunami(%) = '+repr(VR[3])+'\n')
    f.write('VR InSAR LOS(%) = '+repr(VR[4])+'\n')
    f.write('RMS static(%) = '+repr(L2data[0])+'\n')
    f.write('RMS displacement(%) = '+repr(L2data[1])+'\n')
    f.write('RMS velocity(%) = '+repr(L2data[2])+'\n')
    f.write('RMS tsunami(%) = '+repr(L2data[3])+'\n')
    f.write('RMS InSAR LOS(%) = '+repr(L2data[4])+'\n')
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
    row=ifault//nstrike #Row number corresponding to this subfault
    column=ifault-(nstrike*row)
    if nstrike<4 or ndip<4:
        print("ERROR: The fault model is too small for Laplacian regualrization. You need a minimum of 4 rows and 4 columns in the model.")
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
    t1syn=syn.stats.starttime
    t1st=st.stats.starttime
    dt=t1syn-t1st
    if dt>=0: #Synthetic starts before data, pad with zeros
        #How many zeros do I need
        npad=dt/st.stats.delta
        z=zeros(int(npad))
        syn.data=r_[z,syn.data]
        syn.stats.starttime=t1st
    else: #Synthetic starts after the waveform, crop to start fo waveform
        syn.trim(t1st)
    #Now deal with end times
    t2syn=syn.stats.endtime
    t2st=st.stats.endtime
    dt=t2syn-t2st
    if dt>=0: #Synthetic ends after data, crop it
        syn.trim(endtime=t2st)
    else: #Syntetic ends before data, throw an error
        print("ERROR: Synthetic end time is before data end time, recompute longer syntehtics please.")
        return 'Error in GF length'
    return syn
        
def gdims(datafiles,nfaults,decimate):
    '''
    Survey the data files to determine what dimension G will be and return a matrix of zeros 
    with the required dimensions
    '''
    
    from obspy import read
    from numpy import zeros
    from mudpy.green import stdecimate
    
    npts=0
    if decimate==None:
        decimate=1
    for k in range(len(datafiles)):
        e=read(datafiles[k]+'.e')
        n=read(datafiles[k]+'.n')
        u=read(datafiles[k]+'.u')
        if e[0].stats.npts==n[0].stats.npts==u[0].stats.npts:
            npts+=e[0].stats.npts
        else:
            print(str(e[0].stats.npts)+' pts in east component')
            print(str(n[0].stats.npts)+' pts in north component')
            print(str(u[0].stats.npts)+' pts in up component')
            print('ERROR: The 3 components of data are not the same length')
            return 'Error in forming G'
    G=zeros([3*npts,nfaults*2])
    return G 

def gdims_tsun(datafiles,nfaults,decimate):
    '''
    Survey the data files to determine what dimension G will be and return a matrix of zeros 
    with the required dimensions
    '''
    
    from obspy import read
    from numpy import zeros
    from mudpy.green import stdecimate
    
    npts=0
    if decimate==None:
        decimate=1
    for k in range(len(datafiles)):
        tsun=read(datafiles[k])
        npts+=tsun[0].stats.npts
    G=zeros([npts,nfaults*2])
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
        print(sta[k])
        print(gffiles[k])
        out=str(sta[k])+'\t'+repr(round(lon[k],6))+'\t'+repr(round(lat[k],6))+'\t'+str(gffiles[k])+'\n'
        f.write(out)
    f.close()

    
            
def epi2subfault(epicenter,source,vr,tr=0):
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
    
  
def d2epi(epicenter,source):
    '''
    Compute distance from subfault to epicenter
    
    IN:
        epicenter: Epicentral coordinates
        source: Matrix of subfault coordinates
        
    OUT:
        d: distance to epicenter
    '''
    from numpy import tile,sin,cos,deg2rad,sqrt
    #Compute distances from epi to subfault by converting to cartesian
    R=6371
    epicenter=tile(epicenter,(len(source),1))
    xepi=(R-epicenter[:,2])*sin(deg2rad(90-epicenter[:,1]))*cos(deg2rad(epicenter[:,0]))
    yepi=(R-epicenter[:,2])*sin(deg2rad(90-epicenter[:,1]))*sin(deg2rad(epicenter[:,0]))
    zepi=(R-epicenter[:,2])*cos(deg2rad(90-epicenter[:,1]))
    x=(R-source[:,2])*sin(deg2rad(90-source[:,1]))*cos(deg2rad(source[:,0]))
    y=(R-source[:,2])*sin(deg2rad(90-source[:,1]))*sin(deg2rad(source[:,0]))
    z=(R-source[:,2])*cos(deg2rad(90-source[:,1]))
    d=sqrt((xepi-x)**2+(yepi-y)**2+(zepi-z)**2) 
    return d

def get_stats(WG,sol,wd,Ls):
    '''
    Compute basic performance metrics of an inversion
    
    IN:
        WG: Dataweights times GFs WG=W*G
        sol: Solution vector from inversion
        wd: Data weights times data vector, wd=W*d
    OUT:
        L2: ||Gm-d||
        Lm: ||Lm||
    '''
    
    from numpy.linalg import norm
    
    wds=WG.dot(sol)
    L2=norm(wds-wd)
    Lm=norm(Ls.dot(sol))
    return L2,Lm
    
    
def get_VR(home,project_name,GF_list,sol,d,ds,decimate,WG,wd):
    '''
    Compute Variance reduction to the data
    
    IN:
        G: GF matrix
        sol: Solution  vector from inversion
        d: data vector
    OUT:
        VR: Variance reduction (%)
    '''
    
    from numpy import genfromtxt,where,r_,nan
    from obspy import read
    from mudpy.green import stdecimate
    from numpy.linalg import norm
    
    print('... calcualting variance reduction...')
    
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],skip_header=1,dtype='U')
    #Separate into its constituent parts (statics,displacaments, velocities, etc...)
    kstart=0
    kend=0
    
    #Calculate weighted syntehtic data
    wds=WG.dot(sol)
    
    #Statics
    kgf=0
    i=where(GF[:,kgf]==1)[0]
    VRstatic=nan
    L2static=nan
    if len(i)>0:
        for ksta in range(len(i)):
            kend+=3
        #Variance reduction
        res=((d[kstart:kend]-ds[kstart:kend])**2)**0.5
        dnorm=(d[kstart:kend]**2)**0.5 #Yes i know this is dumb, shush
        VRstatic=(1-(res.sum()/dnorm.sum()))*100
        
        #Get RMS (L2 norm) of weigtehd data misfit
        L2static=norm(wds[kstart:kend]-wd[kstart:kend])
   
    #Displacement
    kstart=kend
    kgf=1
    i=where(GF[:,kgf]==1)[0]
    VRdisp=nan
    L2disp=nan
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            n=read(GFfiles[i[ksta],kgf]+'.n')
            if decimate != None:
                n[0]=stdecimate(n[0],decimate)
            npts=n[0].stats.npts
            kend+=3*npts
        #Variance reduction
        res=((d[kstart:kend]-ds[kstart:kend])**2)**0.5
        dnorm=(d[kstart:kend]**2)**0.5 #Yes i know this is dumb, shush
        VRdisp=(1-(res.sum()/dnorm.sum()))*100
        
        #Get RMS (L2 norm) of weigtehd data misfit
        L2disp=norm(wds[kstart:kend]-wd[kstart:kend])
        
    #Velocity
    kstart=kend
    kgf=2
    i=where(GF[:,kgf]==1)[0]
    VRvel=nan
    L2vel=nan
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            n=read(GFfiles[i[ksta],kgf]+'.n')
            if decimate != None:
                n[0]=stdecimate(n[0],decimate)
            npts=n[0].stats.npts
            kend+=3*npts
        #Variance reduction
        res=((d[kstart:kend]-ds[kstart:kend])**2)**0.5
        dnorm=(d[kstart:kend]**2)**0.5 #Yes i know this is dumb, shush
        VRvel=(1-(res.sum()/dnorm.sum()))*100
        
        #Get RMS (L2 norm) of weigtehd data misfit
        L2vel=norm(wds[kstart:kend]-wd[kstart:kend])
        
    #Tsunami
    kstart=kend
    kgf=3
    i=where(GF[:,kgf]==1)[0]
    VRtsun=nan
    L2tsun=nan
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            tsun=read(GFfiles[i[ksta],kgf])
            npts=tsun[0].stats.npts
            kend+=npts
        #Variance reduction
        res=((d[kstart:kend]-ds[kstart:kend])**2)**0.5
        dnorm=(d[kstart:kend]**2)**0.5 #Yes i know this is dumb, shush
        VRtsun=(1-(res.sum()/dnorm.sum()))*100
        
        #Get RMS (L2 norm) of weigtehd data misfit
        L2tsun=norm(wds[kstart:kend]-wd[kstart:kend])
        
    # InSAR
    kstart=kend
    kgf=4
    i=where(GF[:,kgf]==1)[0]
    VRinsar=nan
    L2insar=nan
    if len(i)>0:
        for ksta in range(len(i)):
            kend+=1
        #Variance reduction
        res=((d[kstart:kend]-ds[kstart:kend])**2)**0.5
        dnorm=(d[kstart:kend]**2)**0.5 #Yes i know this is dumb, shush
        VRinsar=(1-(res.sum()/dnorm.sum()))*100  
        
        #Get RMS (L2 norm) of weigtehd data misfit
        L2insar=norm(wds[kstart:kend]-wd[kstart:kend])
        
        
    VR=r_[VRstatic,VRdisp,VRvel,VRtsun,VRinsar]
    L2=r_[L2static,L2disp,L2vel,L2tsun,L2insar]
    
    return VR,L2
    


def get_RMS(home,project_name,run_name,run_number,GF_list,bandpass,use_weights=True):
    '''
    Compute RMS per data type
    
    IN:
        G: GF matrix
        sol: Solution  vector from inversion
        d: data vector
    OUT:
        RMS
    '''
    
    from numpy import genfromtxt,where,r_,nan,zeros,ones
    from obspy import read
    from mudpy.green import stdecimate
    from mudpy.forward import lowpass as lfilter
    
    print('... calcualting RMS...')
    
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    all_weights=genfromtxt(gf_file,usecols=[13,14,15,16,17,18,19,20,21,22,23],skip_header=1,dtype='f8')
    
    # What are the fitler pass bands?
    displacement_bandpass=bandpass[0]
    velocity_bandpass=bandpass[1]
    tsunami_bandpass=bandpass[2]
    
    
    ####    Statics    ####
    kgf=0
    i=where(GF[:,kgf]==1)[0]
    RMSstatic=nan
    datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    
    weights=all_weights[i,0:3]
    nweight=weights[:,0]
    eweight=weights[:,1]
    uweight=weights[:,2]
    
    #where will the data and synthetics go
    n=zeros(len(i))
    e=zeros(len(i))
    u=zeros(len(i))
    ns=zeros(len(i))
    es=zeros(len(i))
    us=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(datapath+str(stations[i[k]])+'.neu')
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        n[k]=neu[0] ; e[k]=neu[1] ; u[k]=neu[2]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.static.neu')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        ns[k]=neus[0] ; es[k]=neus[1] ; us[k]=neus[2]
        
    #Into a single vector they go
    if use_weights==True:
        d_stat=r_[n/nweight,e/eweight,u/uweight]
        ds_stat=r_[ns/nweight,es/eweight,us/uweight]        
    else:
        d_stat=r_[n,e,u]
        ds_stat=r_[ns,es,us]
    RMSstatic=(sum((d_stat-ds_stat)**2)/len(d_stat))**0.5
    
    
    ###  displacememnt waveforms   ####
    kgf=1
    i=where(GF[:,kgf]==1)[0]
    RMSdisplacement=nan
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'  
    datasuffix='disp'
    synthsuffix='disp'
    
    weights=all_weights[i,3:6]
    nweight=weights[:,0]
    eweight=weights[:,1]
    uweight=weights[:,2]
    
    for k in range(len(i)):
        n=read(datapath+stations[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+stations[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+stations[i[k]]+'.'+datasuffix+'.u')
        ns=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.n.sac')
        es=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.e.sac')
        us=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.u.sac')
        
        #Lowpass
        fsample=1./e[0].stats.delta
        e[0].data=lfilter(e[0].data,displacement_bandpass,fsample,2)
        n[0].data=lfilter(n[0].data,displacement_bandpass,fsample,2)
        u[0].data=lfilter(u[0].data,displacement_bandpass,fsample,2)
        es[0].data=lfilter(es[0].data,displacement_bandpass,fsample,2)
        ns[0].data=lfilter(ns[0].data,displacement_bandpass,fsample,2)
        us[0].data=lfilter(us[0].data,displacement_bandpass,fsample,2)
        
        if use_weights==True:
            nw=ones(n[0].stats.npts)*nweight[k]
            ew=ones(n[0].stats.npts)*eweight[k]
            uw=ones(n[0].stats.npts)*uweight[k]
            dtemp=r_[n[0].data/nw,e[0].data/ew,u[0].data/uw]
            dstemp=r_[ns[0].data/nw,es[0].data/ew,us[0].data/uw]  
        else:
            dtemp=r_[n[0].data,e[0].data,u[0].data]
            dstemp=r_[ns[0].data,es[0].data,us[0].data]
        
        if k==0:
            d_disp=dtemp.copy()
            ds_disp=dstemp.copy()
        else:
            d_disp=r_[d_disp,dtemp.copy()]
            ds_disp=r_[ds_disp,dstemp.copy()]
    
    RMSdisplacement=(sum((d_disp-ds_disp)**2)/len(d_disp))**0.5
    
    
    ###  velocity waveforms   ####
    kgf=2
    i=where(GF[:,kgf]==1)[0]
    RMSvelocity=nan
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'  
    datasuffix='vel'
    synthsuffix='vel'
    
    weights=all_weights[i,6:9]
    nweight=weights[:,0]
    eweight=weights[:,1]
    uweight=weights[:,2]
    
    for k in range(len(i)):
        n=read(datapath+stations[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+stations[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+stations[i[k]]+'.'+datasuffix+'.u')
        ns=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.n.sac')
        es=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.e.sac')
        us=read(synthpath+run_name+'.'+run_number+'.'+stations[i[k]]+'.'+synthsuffix+'.u.sac')
        
        #Bandpass
        fsample=1./e[0].stats.delta
        e[0].data=lfilter(e[0].data,velocity_bandpass,fsample,2)
        n[0].data=lfilter(n[0].data,velocity_bandpass,fsample,2)
        u[0].data=lfilter(u[0].data,velocity_bandpass,fsample,2)
        es[0].data=lfilter(es[0].data,velocity_bandpass,fsample,2)
        ns[0].data=lfilter(ns[0].data,velocity_bandpass,fsample,2)
        us[0].data=lfilter(us[0].data,velocity_bandpass,fsample,2)
        
        if use_weights==True:
            nw=ones(n[0].stats.npts)*nweight[k]
            ew=ones(n[0].stats.npts)*eweight[k]
            uw=ones(n[0].stats.npts)*uweight[k]
            dtemp=r_[n[0].data/nw,e[0].data/ew,u[0].data/uw]
            dstemp=r_[ns[0].data/nw,es[0].data/ew,us[0].data/uw]  
        else:
            dtemp=r_[n[0].data,e[0].data,u[0].data]
            dstemp=r_[ns[0].data,es[0].data,us[0].data]
        
        if k==0:
            d_vel=dtemp.copy()
            ds_vel=dstemp.copy()
        else:
            d_vel=r_[d_vel,dtemp.copy()]
            ds_vel=r_[ds_vel,dstemp.copy()]
    
    RMSvelocity=(sum((d_vel-ds_vel)**2)/len(d_vel))**0.5
       
    RMStsunami=nan
    RMSinsar=nan
    
    # Total RMS
    d=r_[d_stat,d_disp,d_vel]
    ds=r_[ds_stat,ds_disp,ds_vel]
    RMStotal=(sum((d-ds)**2)/len(d))**0.5
    
    RMS=r_[RMSstatic,RMSdisplacement,RMSvelocity,RMStsunami,RMSinsar,RMStotal]
    
    return RMS


    
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
        print("... Static ABIC requested")
        s=norm(d-G.dot(sol))**2+(lambda_s**2)*norm(Ls.dot(sol))**2
        a1=N*log(s)
        a2=M*log(lambda_s**2)
        sq,a3=slogdet(GTG+(lambda_s**2)*LsLs)
        #Add 'em up
        ABIC=a1-a2+a3
        return ABIC
    else: #There is a double regularization, use Fukahata et al. definition
        print('... computing 2d-ABIC')
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
    from mudpy.forward import get_mu
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Open structure file
    mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
    #Get slip quantities
    iss=2*arange(len(sol)//2)
    ids=2*arange(len(sol)//2)+1
    ss=sol[iss]
    ds=sol[ids]
    #Get total slip
    slip=(ss**2+ds**2)**0.5
    #Get rigidities and areas
    mu=zeros(len(ds))
    A=zeros(len(ds))
    i=0
    for krupt in range((len(sol)//len(f))//2):
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
    from numpy import array,deg2rad,cos,sin,arange,hstack,zeros,expand_dims
    
    #Split into strike-slip and dip-slip
    iss=2*arange(0,len(sol)/2,1)
    ids=2*arange(0,len(sol)/2,1)+1
    iss=iss.astype('int')
    ids=ids.astype('int')
    if len(iss)==1:
        ss=sol[0]
        ds=sol[1]
    else:
        ss=sol[iss]
        ds=sol[ids]
    #Rotate
    beta=deg2rad(beta)
    rot=array([[cos(beta),sin(beta)],[-sin(beta),cos(beta)]]).dot(hstack((ss,ds)).T)
    #Re-insert in output vector
    out=zeros(sol.shape)
    out[iss]=expand_dims(rot[0,:].T,1)
    out[ids]=expand_dims(rot[1,:].T,1)
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
    from numpy import interp
    from mudpy.forward import round_time
    
    #Get data sampling interval
    delta=data.stats.delta
    #Get synthetic start time
    t1synth=synth.stats.starttime
    #Get start time sampled at data interval and correction encessary
    t1=round_time(t1synth,delta)
    dt=t1-t1synth #This is the correctiont hat needs to be applied
    #Make current time vector
    tsynth=synth.times()
    tcorrect=tsynth+dt
    #Interpolate to new time vector
    synth_data=synth.data
    synth_correct=interp(tcorrect,tsynth,synth_data)
    #Place in output stream object
    synth.starttime=t1
    synth.data=synth_correct
    return synth

def resample_synth_tsun(synth,data):
    '''
    Resample a synthetic to the time samples contained in the input data
    IN:
        synth: synthetic stream object
        data: data stream object
    OUT:
        st: resampled synthetic
    '''
    from numpy import interp,arange
    from mudpy.forward import round_time
    
    #Get data sampling interval
    delta=data.stats.delta
    #Get synthetic start time
    t1synth=synth.stats.starttime
    #Get start time sampled at data interval and correction encessary
    t1=round_time(t1synth,delta)
    dt=t1-t1synth #This is the correctiont hat needs to be applied
    #Make current time vector
    tsynth=synth.times()
    tcorrect=arange(tsynth[0],tsynth[-1],delta)+dt
    #Interpolate to new time vector
    synth_data=synth.data
    synth_correct=interp(tcorrect,tsynth,synth_data)
    #Place in output stream object
    synth.starttime=t1
    synth.data=synth_correct
    synth.stats.delta=delta
    return synth
    
def model_covariance(home,project_name,run_name,run_number,fault_name,G_name,nfaults,
                        num_windows,bounds,GF_list,decimate,lowpass,beta):
    '''
    Compute model covariance matrix
    '''
    from numpy import load,arange,zeros,genfromtxt,expand_dims,save
    from numpy.linalg import norm,inv
    from datetime import datetime
    import gc

    t0=datetime.now()
    #Get smoothing from log file
    outdir=home+project_name+'/output/inverse_models/models/'
    log=outdir+run_name+'.'+run_number+'.log'
    gf_file=home+project_name+'/data/station_info/'+GF_list
    #Get value of smoothign parameter
    with open(log) as f:
        for line in f:
            if 'lambda_spatial' in line:
                ls=float(line.split('=')[1])
            if 'lambda_temporal' in line:
                lt=float(line.split('=')[1])
    #Read data variances and make data covariance
    #Dos trike slip and dip slip in different stages
    print('Getting data covariance...')
    #Cd=data_covariance(gf_file,decimate)
    #Cd=diag(Cd)
    #Cd=diag(1/Cd)
    #Load G
    print('Computing for SS model parameters')
    print('Getting G...')
    G_name=home+project_name+'/GFs/matrices/'+G_name
    print('Loading '+G_name)
    if G_name[-3:]!='npy':
            G_name=G_name+'.npy'
    G=load(G_name)
    iss=2*arange(G.shape[1]/2)
    ids=2*arange(G.shape[1]/2)+1
    gc.collect()
    print('Compressing G...')
    #Get regularization matrices
    print('Getting regularization matrices...')
    Ls=getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
    Lt=getLt(home,project_name,fault_name,num_windows)
    #Read model
    print('Reading model')
    model=outdir+run_name+'.'+run_number+'.inv'
    model=genfromtxt(model,usecols=(8,9))
    m0=zeros((model.size,1))
    m0[iss]=expand_dims(model[:,0],1) 
    m0[ids]=expand_dims(model[:,1],1) 
    #rotate
    m0=ds2rot(m0,beta)
    #Get data
    print('Reading data...')
    d=getdata(home,project_name,GF_list,decimate,lowpass)
    #Compute residual
    sm0=norm(d-G.dot(m0))**2+(ls**2)*norm(Ls.dot(m0))**2+(lt**2)*norm(Lt.dot(m0))**2
    #Get sigma
    sigma=sm0/len(d)
    #Compute model covariance
    print('Computing model covariance...')
    #Cd=csc_matrix(Cd)
    Cm=(G.T).dot(G)+(ls**2)*(Ls.T).dot(Ls)+(lt**2)*(Lt.T).dot(Lt)
    print('Inverting G\'G')
    Cm=inv(Cm)
    Cm=Cm*sigma
    print('Extracting diagonal terms...')
    Cm=Cm.diagonal()
    G=None
    #Prep for output
    Cm_name=G_name.split('.npy')[0]+'.cov'
    print('Saving to '+Cm_name)
    save(Cm_name,Cm)
    deltaT=datetime.now()-t0
    print('Well that only took '+str(deltaT))
    
    
    
def data_covariance(gf_file,decimate):
    '''
    Form data covariance matrix
    '''
    from numpy import where,genfromtxt,r_,diag,ones
    from obspy import read
    from mudpy.green import stdecimate
    
    gflist=genfromtxt(gf_file,usecols=(3,4,5,6,7))
    data_files=genfromtxt(gf_file,usecols=(9,10),dtype='S')
    #Get static covariance
    sigma_static=[]
    i=where(gflist[:,0]==1)[0]
    if len(i)>0:
        sigman=genfromtxt(gf_file,usecols=13)
        sigman=sigman[i]
        sigmae=genfromtxt(gf_file,usecols=14)
        sigmae=sigmae[i]
        sigmau=genfromtxt(gf_file,usecols=15)
        sigmau=sigmau[i]
        sigma_static=r_[sigman,sigmae,sigmau]
    #Get displacement covariance
    sigma_disp=[]
    i=where(gflist[:,1]==1)[0]
    if len(i)>0:
        #Read variances
        sn=genfromtxt(gf_file,usecols=16)
        se=genfromtxt(gf_file,usecols=17)
        su=genfromtxt(gf_file,usecols=18)
        for k in range(len(i)):
            #Get length of time series
            st=read(data_files[i[k],0]+'.n')
            st[0]=stdecimate(st[0],decimate)
            Nt=st[0].stats.npts
            #Cocnatenate
            if k==0:
                sigman=ones(Nt)*sn[i[k]]
                sigmae=ones(Nt)*se[i[k]]
                sigmau=ones(Nt)*su[i[k]]
            else:
                sigman=r_[sigman,ones(Nt)*sn[i[k]]]
                sigmae=r_[sigmae,ones(Nt)*se[i[k]]]
                sigmau=r_[sigmau,ones(Nt)*su[i[k]]]
        sigma_disp=r_[sigman,sigmae,sigmau]
    #Get velocioty covariance
    sigma_vel=[]
    i=where(gflist[:,2]==1)[0]
    if len(i)>0:
        #Read variances
        sn=genfromtxt(gf_file,usecols=19)
        se=genfromtxt(gf_file,usecols=20)
        su=genfromtxt(gf_file,usecols=21)
        for k in range(len(i)):
            #Get length of time series
            st=read(data_files[i[k],1]+'.n')
            st[0]=stdecimate(st[0],decimate)
            Nt=st[0].stats.npts
            #Cocnatenate
            if k==0:
                sigman=ones(Nt)*sn[i[k]]
                sigmae=ones(Nt)*se[i[k]]
                sigmau=ones(Nt)*su[i[k]]
            else:
                sigman=r_[sigman,ones(Nt)*sn[i[k]]]
                sigmae=r_[sigmae,ones(Nt)*se[i[k]]]
                sigmau=r_[sigmau,ones(Nt)*su[i[k]]]
        sigma_vel=r_[sigman,sigmae,sigmau]
    Cd=diag(r_[sigma_static,sigma_disp,sigma_vel])
    return Cd
    
    
    
def make_tgf_dtopo(home,project_name,model_name,tgf_file,coast_file,fault_name,
            time_epi,tsun_dt,maxt,topo_dx_file,topo_dy_file,
            topo_effect=False,instantaneous=True,hot_start=0,average_instantaneous=10):
    '''
    Create moving topography input files for geoclaw
    
    tsun_dt - Sampling itnerval of synthetics ALREADY made
    model_name is structure file
    time_epi is UTC string
    
    '''
    import datetime
    from numpy import genfromtxt,zeros,arange,meshgrid,ones,c_,savetxt,argmin,nan,where,mean
    from obspy import read
    from netCDF4 import Dataset
    from scipy.interpolate import griddata

    #Get station names
    staname=genfromtxt(home+project_name+'/data/station_info/'+tgf_file,usecols=0,dtype='U')
    sta=genfromtxt(home+project_name+'/data/station_info/'+tgf_file)
    lon=sta[:,1]
    lat=sta[:,2]
    #Get fault file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Where is the data
    green_dir=home+project_name+'/GFs/tsunami/'
    #green_dir=home+project_name+'/GFs/dynamic/'
    #Define time deltas
    td_max=datetime.timedelta(seconds=maxt)
    td=datetime.timedelta(seconds=tsun_dt)
    #Maximum time to be modeled
    tmax=time_epi+td_max
    Nt=(tmax-time_epi)/tsun_dt
    #Read topo/bathy derivative fileksub=0
    #if topo_dx_file!=None
    bathy_dx=Dataset(topo_dx_file,'r')
    zdx=bathy_dx.variables['z'][:]
    bathy_dy=Dataset(topo_dy_file,'r')
    zdy=bathy_dy.variables['z'][:]
    #Make mesh that matches it for interpolation later
    try:
        lonb=bathy_dx.variables['lon'][:]
        latb=bathy_dx.variables['lat'][:]
    except:
        lonb=bathy_dx.variables['x'][:]
        latb=bathy_dx.variables['y'][:]
    #print("Correcting longitude"
    #lon=360+lon
    #lonb=360+lonb
    loni,lati=meshgrid(lonb,latb)
    #Now apply motion from every subfault
    for ksub in range(hot_start,len(f)): #Loops through subfaults
        if ksub%10==0:
            print('... working on subfault '+str(ksub)+' of '+str(len(f)))
        #Depth string
        zs=f[ksub,3]
        strdepth='%.4f' % zs
        #subfault number
        sub=str(ksub+1).rjust(4,'0')
        #Subfault dir
        subdir=green_dir+model_name+'_'+strdepth+'.sub'+sub+'/'
        for ksta in range(len(sta)):
            if ksta%500==0:
                print('... ... working on seafloor grid point '+str(ksta)+' of '+str(len(sta)))
            eds=read(subdir+staname[ksta]+'.subfault'+sub+'.DS.disp.e')
            nds=read(subdir+staname[ksta]+'.subfault'+sub+'.DS.disp.n')
            uds=read(subdir+staname[ksta]+'.subfault'+sub+'.DS.disp.z')
            ess=read(subdir+staname[ksta]+'.subfault'+sub+'.SS.disp.e')
            nss=read(subdir+staname[ksta]+'.subfault'+sub+'.SS.disp.n')
            uss=read(subdir+staname[ksta]+'.subfault'+sub+'.SS.disp.z')
            eds=interp_and_resample(eds,1.0,time_epi)
            nds=interp_and_resample(nds,1.0,time_epi)
            uds=interp_and_resample(uds,1.0,time_epi)
            ess=interp_and_resample(ess,1.0,time_epi)
            nss=interp_and_resample(nss,1.0,time_epi)
            uss=interp_and_resample(uss,1.0,time_epi)
            #Keep only data between time_epi and tmax
            eds.trim(time_epi,tmax,fill_value=eds[0].data[-1],pad=True)
            nds.trim(time_epi,tmax,fill_value=nds[0].data[-1],pad=True)
            uds.trim(time_epi,tmax,fill_value=uds[0].data[-1],pad=True)
            ess.trim(time_epi,tmax,fill_value=ess[0].data[-1],pad=True)
            nss.trim(time_epi,tmax,fill_value=nss[0].data[-1],pad=True)
            uss.trim(time_epi,tmax,fill_value=uss[0].data[-1],pad=True)
            #Decimate to original smapling interval
            eds[0].decimate(2,no_filter=True)
            nds[0].decimate(2,no_filter=True)
            uds[0].decimate(2,no_filter=True)
            ess[0].decimate(2,no_filter=True)
            nss[0].decimate(2,no_filter=True)
            uss[0].decimate(2,no_filter=True)
            #Initalize matrices
            if ksta==0:
                eDS=zeros((nss[0].stats.npts,len(sta)))
                eSS=eDS.copy()
                nDS=eDS.copy()
                nSS=eDS.copy()
                uDS=eDS.copy()
                uSS=eDS.copy()
            #Populate matrix
            eDS[:,ksta]=eds[0].data
            nDS[:,ksta]=nds[0].data
            uDS[:,ksta]=uds[0].data
            eSS[:,ksta]=ess[0].data
            nSS[:,ksta]=nss[0].data
            uSS[:,ksta]=uss[0].data
        #Only apply topo effect to some points
        if coast_file!=None:
            #Straight line coordinates
            coast=genfromtxt(coast_file)
            #Get mask for applying horizontal effect
            mask=zeros(loni.shape)
            for k1 in range(loni.shape[0]):
                for k2 in range(loni.shape[1]):
                    #if (lati[k1,k2]-b)/m>loni[k1,k2]: #Point is to the left (or right), do not apply horizontal effect
                    #Find two closest lat,lon points
                    ip=argmin(abs(lati[k1,k2]-coast[:,1]))
                    if loni[k1,k2]>coast[ip,0]:
                        mask[k1,k2]=nan
            imask1,imask2=where(mask==0)#Points tot he right DO apply horiz. effect
        #Now go one epoch at a time, and interpolate all fields
        print('... interpolating coseismic offsets to a regular grid')
        if instantaneous==True:
            nt_iter=2
        else:
            nt_iter=uSS.shape[0]
        for kt in range(nt_iter):
            if kt%20==0:
                print('... ... working on time slice '+str(kt)+' of '+str(nt_iter))
            if instantaneous==True and kt==1:
                
                nds=griddata((lon,lat),mean(nDS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
                eds=griddata((lon,lat),mean(eDS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
                uds=griddata((lon,lat),mean(uDS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
                nss=griddata((lon,lat),mean(nSS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
                ess=griddata((lon,lat),mean(eSS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
                uss=griddata((lon,lat),mean(uSS[-average_instantaneous:,:],axis=0),(loni,lati),method='linear',fill_value=0)
            else:
                nds=griddata((lon,lat),nDS[kt,:],(loni,lati),method='linear',fill_value=0)
                eds=griddata((lon,lat),eDS[kt,:],(loni,lati),method='linear',fill_value=0)
                uds=griddata((lon,lat),uDS[kt,:],(loni,lati),method='linear',fill_value=0)
                nss=griddata((lon,lat),nSS[kt,:],(loni,lati),method='linear',fill_value=0)
                ess=griddata((lon,lat),eSS[kt,:],(loni,lati),method='linear',fill_value=0)
                uss=griddata((lon,lat),uSS[kt,:],(loni,lati),method='linear',fill_value=0)
            #Output vertical
            uout_ds=uds
            uout_ss=uss
            #Apply effect of topography advection
            if topo_effect==True:
                uout_ds[imask1,imask2]=uout_ds[imask1,imask2]+zdx[imask1,imask2]*eds[imask1,imask2]+zdy[imask1,imask2]*nds[imask1,imask2]
                uout_ss[imask1,imask2]=uout_ss[imask1,imask2]+zdx[imask1,imask2]*ess[imask1,imask2]+zdy[imask1,imask2]*nss[imask1,imask2]
            #Convert to column format and append
            xyz_ds=grd2xyz(uout_ds,loni,lati)
            xyz_ss=grd2xyz(uout_ss,loni,lati)
            tvec=(kt*tsun_dt)*ones((len(xyz_ds),1))
            if kt==0: #Intialize
                numel=uout_ds.size #Number of elements in grid
                kwrite=numel #Where to write the data
                dtopo_ds=zeros((numel*nt_iter,4))
                dtopo_ss=zeros((numel*nt_iter,4))
                dtopo_ds[0:kwrite,1:3]=xyz_ds[:,0:2]
                dtopo_ss[0:kwrite,1:3]=xyz_ss[:,0:2]
            else:
                if instantaneous==True:
                    tvec=ones(len(tvec))
                dtopo_ds[kwrite:kwrite+numel,:]=c_[tvec,xyz_ds]
                dtopo_ss[kwrite:kwrite+numel,:]=c_[tvec,xyz_ss]
                kwrite=kwrite+numel
        print('... writting dtopo files')
        savetxt(subdir+'DS.dtopo',dtopo_ds,fmt='%i\t%.6f\t%.6f\t%.4e')
        savetxt(subdir+'SS.dtopo',dtopo_ss,fmt='%i\t%.6f\t%.6f\t%.4e')
             
                
def make_tgf_dtopo_static(home,project_name,model_name,tgf_file,fault_name,
            time_epi,tsun_dt,maxt,topo_dx_file,topo_dy_file,
            topo_effect=False,hot_start=0):
    '''
    Create moving topography input files for geoclaw
    
    tsun_dt - Sampling itnerval of synthetics ALREADY made
    model_name is structure file
    time_epi is UTC string
    
    '''
    import datetime
    from numpy import genfromtxt,zeros,arange,meshgrid,ones,c_,savetxt,argmin,nan,where,mean,expand_dims
    from obspy import read
    from netCDF4 import Dataset
    from scipy.interpolate import griddata

    #Get station names
    staname=genfromtxt(home+project_name+'/data/station_info/'+tgf_file,usecols=0,dtype='U')
    sta=genfromtxt(home+project_name+'/data/station_info/'+tgf_file)
    lon=sta[:,1]
    lat=sta[:,2]
    #Get fault file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Where is the data
    green_dir=home+project_name+'/GFs/tsunami/'

    #Read topo/bathy derivative fileksub=0
    #if topo_dx_file!=None
    bathy_dx=Dataset(topo_dx_file,'r')
    zdx=bathy_dx.variables['z'][:]
    bathy_dy=Dataset(topo_dy_file,'r')
    zdy=bathy_dy.variables['z'][:]
    #Make mesh that matches it for interpolation later
    try:
        lonb=bathy_dx.variables['lon'][:]
        latb=bathy_dx.variables['lat'][:]
    except:
        lonb=bathy_dx.variables['x'][:]
        latb=bathy_dx.variables['y'][:]
    #print("Correcting longitude"
    #lon=360+lon
    #lonb=360+lonb
    loni,lati=meshgrid(lonb,latb)
    #Now apply motion from every subfault
    for ksub in range(hot_start,len(f)): #Loops through subfaults
        if ksub%10==0:
            print('... working on subfault '+str(ksub)+' of '+str(len(f)))
        #Depth string
        zs=f[ksub,3]
        strdepth='%.4f' % zs
        #subfault number
        sub=str(ksub+1).rjust(4,'0')
        #Subfault dir
        subdir=green_dir+model_name+'_'+strdepth+'.sub'+sub+'/'
            
        ds=genfromtxt(subdir+'subfault'+sub+'.DS.static.neu')
        ss=genfromtxt(subdir+'subfault'+sub+'.SS.static.neu')
        
        eds=ds[:,2]
        nds=ds[:,1]
        uds=ds[:,3]
        ess=ss[:,2]
        nss=ss[:,1]
        uss=ss[:,3]


        #Populate matrix
        eDS=expand_dims(eds,1).T
        nDS=expand_dims(nds,1).T
        uDS=expand_dims(uds,1).T
        eSS=expand_dims(ess,1).T
        nSS=expand_dims(nss,1).T
        uSS=expand_dims(uss,1).T

        #Now go one epoch at a time, and interpolate all fields
        print('... ... interpolating coseismic offsets to a regular grid')
        
        
        nds=griddata((lon,lat),nDS[0,:],(loni,lati),method='linear',fill_value=0)
        eds=griddata((lon,lat),eDS[0,:],(loni,lati),method='linear',fill_value=0)
        uds=griddata((lon,lat),uDS[0,:],(loni,lati),method='linear',fill_value=0)
        nss=griddata((lon,lat),nSS[0,:],(loni,lati),method='linear',fill_value=0)
        ess=griddata((lon,lat),eSS[0,:],(loni,lati),method='linear',fill_value=0)
        uss=griddata((lon,lat),uSS[0,:],(loni,lati),method='linear',fill_value=0)
            
            
        #Output vertical
        uout_ds=uds
        uout_ss=uss
        
 
        #Convert to column format and append
        xyz_ds=grd2xyz(uout_ds,loni,lati)
        xyz_ss=grd2xyz(uout_ss,loni,lati)

        numel=uout_ds.size #Number of elements in grid
        dtopo_ds=zeros((numel*2,4))
        dtopo_ss=zeros((numel*2,4))
        dtopo_ds[0:numel,1:3]=xyz_ds[:,0:2]
        dtopo_ss[0:numel,1:3]=xyz_ss[:,0:2]
            
        dtopo_ds[numel:,0]=1 
        dtopo_ss[numel:,0]=1            
        dtopo_ds[numel:,1:4]=xyz_ds
        dtopo_ss[numel:,1:4]=xyz_ss

        print('... ... writting dtopo files')
        savetxt(subdir+'DS.dtopo',dtopo_ds,fmt='%i\t%.6f\t%.6f\t%.4e')
        savetxt(subdir+'SS.dtopo',dtopo_ss,fmt='%i\t%.6f\t%.6f\t%.4e')
        
        
        
                      
def tsunami_gf(home,project_name,model_name,fault_name,hot_start):
    '''
    Create Geoclaw parameter files and make system call to run simulation
    '''
    from numpy import genfromtxt,r_,arange
    from os import makedirs,chdir
    from shutil import copy
    from shlex import split
    import subprocess
    import gc

    #Load fault file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Where is the data
    green_dir=home+project_name+'/GFs/tsunami/'  
    for ksub in range(hot_start,len(f)):
        print('... working on subfault '+str(ksub)+' of '+str(len(f)))
        
        #Depth string
        zs=f[ksub,3]        
        strdepth='%.4f' % zs
        #subfault number
        sub=str(int(f[ksub,0])).rjust(4,'0')
        #Subfault dir
        subdir=green_dir+model_name+'_'+strdepth+'.sub'+sub+'/'
        #Make dirs, one for DS one for SS
        try:
            makedirs(subdir+'_DS')
        except:
            print('Warning: DS directory already exists')
        try:
            makedirs(subdir+'_SS')
        except:
            print('Warning: SS directory already exists')
        #Move setrun.py and makefile into each directory
        copy(home+project_name+'/data/station_info/setrun.py',subdir+'_DS/')
        copy(home+project_name+'/data/station_info/Makefile',subdir+'_DS/')
        copy(home+project_name+'/data/station_info/setrun.py',subdir+'_SS/')
        copy(home+project_name+'/data/station_info/Makefile',subdir+'_SS/')
        #Modify the dtopo declaration
        dtopoDS=subdir+'DS.dtopo'
        dtopoSS=subdir+'SS.dtopo'
        s_ds = open(subdir+'_DS/setrun.py').read()
        s_ss = open(subdir+'_SS/setrun.py').read()
        s_ds = s_ds.replace('GF_dtopofile',dtopoDS)
        s_ss = s_ss.replace('GF_dtopofile',dtopoSS)  
        fpy = open(subdir+'_DS/setrun.py', 'w')
        fpy.write(s_ds)
        fpy.close()
        fpy = open(subdir+'_SS/setrun.py', 'w')
        fpy.write(s_ss)
        fpy.close()
        #Now run it
        geoclawDS='make .output'
        geoclawSS='make .output'
        chdir(subdir+'_DS')
        subprocess.call(split(geoclawDS))
        #p.communicate()
        gc.collect()
        chdir(subdir+'_SS')
        subprocess.call(split(geoclawSS))
        #p.communicate() 
        gc.collect()
        
                  
def tsunami2sac(home,project_name,model_name,fault_name,tlims,dt,time_epi,hot_start):
    '''
    Create Geoclaw parameter files and make system call to run simulation
    '''
    from numpy import genfromtxt,where,intersect1d,arange,interp,r_
    from obspy import Stream,Trace
    
    #Load fault file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Load gauge correspondences
    gaugesGC=genfromtxt(home+project_name+'/data/station_info/gauges.dict',usecols=0,dtype='U')  #This is what they are called in GeoClaw
    gauges=genfromtxt(home+project_name+'/data/station_info/gauges.dict',usecols=3,dtype='U')   #their actual name
    #Where is the data
    green_dir=home+project_name+'/GFs/tsunami/'  
    for ksub in range(hot_start,len(f)):
        #Depth string
        zs=f[ksub,3]        
        strdepth='%.4f' % zs
        #subfault number
        sub=str(int(f[ksub,0])).rjust(4,'0')
        #Subfault dir
        subdir=green_dir+model_name+'_'+strdepth+'.sub'+sub+'/'
        print('... working on  '+subdir)
        for kgauge in range(gauges.size):
            #Read gauge.fort file
            if gauges.size<2:
                dsfile=subdir+'_DS/_output/gauge%s.txt' % (gaugesGC[kgauge])
                ssfile=subdir+'_SS/_output/gauge%s.txt' % (gaugesGC[kgauge])
            else:
                dsfile=subdir+'_DS/_output/gauge%s.txt' % (gaugesGC[kgauge])
                ssfile=subdir+'_SS/_output/gauge%s.txt' % (gaugesGC[kgauge])
            gds=genfromtxt(dsfile)
            gss=genfromtxt(ssfile)
            
            #Get time vector
            tss=gss[:,1]
            tds=gds[:,1]
            #Get data
            etass=gss[:,5]
            etads=gds[:,5]
            #Trim and resample to regular interval
            itrim=intersect1d(where(tss>=tlims[0])[0],where(tss<=tlims[1])[0])
            tss=tss[itrim]
            etass=etass[itrim]
            itrim=intersect1d(where(tds>=tlims[0])[0],where(tds<=tlims[1])[0])
            tds=tds[itrim]
            etads=etads[itrim]
            ti=arange(tlims[0],tlims[1]+dt,dt)
            etass=interp(ti,tss,etass)
            etads=interp(ti,tds,etads)
            tiout=arange(0,tlims[1]+dt,dt)
            etass=interp(tiout,r_[0,ti],r_[0,etass])
            etads=interp(tiout,r_[0,ti],r_[0,etads])
            #etass=lowpass(etass,1./dt,fcorner,4)

            
            #Initalize stream objects
            stds=Stream(Trace())
            stss=Stream(Trace())
            #Populate fields
            stss[0].data=etass
            stss[0].stats.delta=dt
            stss[0].stats.starttime=time_epi
            stds[0].data=etads
            stds[0].stats.delta=dt
            stds[0].stats.starttime=time_epi
            #Write to file
            if gauges.size<2:
                stss.write(subdir+str(gauges)+'.ss.tsun',format='SAC')
                stds.write(subdir+str(gauges)+'.ds.tsun',format='SAC')
            else:
                stss.write(subdir+str(gauges[kgauge])+'.ss.tsun',format='SAC')
                stds.write(subdir+str(gauges[kgauge])+'.ds.tsun',format='SAC')


    
def grd2xyz(uout,lon,lat):
    ''''
    Write dtopo file
    '''
    from numpy import zeros
    nlat=uout.shape[0]
    nlon=uout.shape[1]
    k=0
    xyz=zeros((nlon*nlat,3))
    for klat in range(1,nlat+1):
        for klon in range(nlon):
            xyz[k,0]=lon[-klat,klon]
            xyz[k,1]=lat[-klat,klon]
            xyz[k,2]=uout[-klat,klon]
            k+=1
    return xyz
    
            
def interp_and_resample(st,dt,time_epi):
    '''
    First interpolate to dt sampling rate, then resmaple so a sample exists at the
    epicentral time
    '''
    from numpy import interp,arange
    from datetime import timedelta
    stout=st.copy()
    t=st[0].times()
    t1=st[0].stats.starttime
    mu=t1.microsecond
    ti=arange(t[0]-mu/1e6,t[-1]-mu/1e6,dt)
    y=st[0].data
    yi=interp(ti,t,y)
    stout[0].data=yi
    stout[0].stats.delta=dt
    stout[0].stats.starttime=t1-timedelta(microseconds=mu)
    return stout
    
    
def data_norms(home,project_name,GF_list,decimate=None,bandpass=[None,None,None]):
    '''
    Read each data type and extract it's norm, this is used to inform the
    weighting scheme
    '''
    from numpy import genfromtxt,where
    from scipy.linalg import norm
    from mudpy.inverse import getdata
    from obspy import read
    
    #Read data vector
    d=getdata(home,project_name,GF_list,decimate,bandpass=bandpass)
    #Figure out which indices belong to which data type
    gf_file=home+project_name+'/data/station_info/'+GF_list
    #Read station flags
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    # how many of each type?
    istatic=where(GF[:,0]==1)[0]
    idisp=where(GF[:,1]==1)[0]
    ivel=where(GF[:,2]==1)[0]
    itsun=where(GF[:,3]==1)[0]
    iinsar=where(GF[:,4]==1)[0]
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='U') 
    kstart=0
    kend=0
    #compute norms
    if len(istatic)>0:
        kend+=len(istatic)*3
        N=norm(d[0:kend])
        print("||statics|| = "+str(N))
    else:
        print("||statics|| = NaN")
    if len(idisp)>0: #read displacememnt waveforms
        kstart=kend
        for kfile in range(len(idisp)):
            n=read(GFfiles[idisp[kfile],1]+'.n')
            kend+=3*n[0].stats.npts
        N=norm(d[kstart:kend])
        print("||disp. waveforms|| = "+str(N))
    else:
        print("||disp. waveforms|| = NaN")
    if len(ivel)>0:
        kstart=kend
        for kfile in range(len(ivel)):
            n=read(GFfiles[ivel[kfile],2]+'.n')
            kend+=3*n[0].stats.npts
        N=norm(d[kstart:kend])
        print("||vel. waveforms|| = "+str(N))
    else:
        print("||vel. waveforms|| = NaN")
    if len(itsun)>0:
        kstart=kend
        for kfile in range(len(itsun)):
            tsun=read(GFfiles[itsun[kfile],3])
            kend+=tsun[0].stats.npts
        N=norm(d[kstart:kend])
        print("||tsunami|| = "+str(N))
    else:
        print("||tsunami|| = NaN")
    if len(iinsar)>0:
        kstart=kend
        kend+=len(iinsar)
        N=norm(d[kstart:kend])
        print("||InSAR|| = "+str(N))
    else:
        print("||InSAR|| = NaN")
    