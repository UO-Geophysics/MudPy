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
        #Done, now concatenate them all and save to binary file
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
            n=dtemp[0]
            e=dtemp[1]
            u=dtemp[2]
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
    
    from numpy import loadtxt,zeros
    
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
            if type(stencil)!=bool: #No errors were reported
                #Add strike slip branches of stencil
                L[2*kfault,2*stencil]=1
                #Add dip slip branches of stencil
                L[2*kfault+1,2*stencil+1]=1
                #Add strike slip central node with correction
                #correction=0
                L[2*kfault,2*kfault]=-4+correction
                #Add dip slip central node with correction
                #correction=0
                L[2*kfault+1,2*kfault+1]=-4+correction
            else:
                return False
        h=zeros(len(L))
        return L,h
    else:
        print 'ERROR: Unknown regularization type ('+regularization_type+') requested.'
        return False
        
def get_data_weights(home,project_name,GF_list,d,rupture_speeds):
    '''
    Assemble matrix of data weights from sigmas of observations
    '''    
    from numpy import genfromtxt,where,zeros,ones,diag_indices_from
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
    for krup in range(len(rupture_speeds)):
        for ksta in range(len(i)):
            #Read waveform to determine length of insert
            st=read(GFfiles[i[ksta],kgf]+'.n')
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
    for krup in range(len(rupture_speeds)):
        for ksta in range(len(i)):
            #Read waveform to determine length of insert
            st=read(GFfiles[i[ksta],kgf]+'.n')
            nsamples=st[0].stats.npts
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
                    ###### These need to be changed to neu, enu is stupid Diego ########
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
    
#=================        Write inversion results      =========================
    
def write_model(home,project_name,run_name,fault_name,model_name,rupture_speeds,epicenter,sol,num):
    '''
    Write model results
    '''
    
    from numpy import genfromtxt,loadtxt,arange,zeros,c_,savetxt
    from forward import get_mu
    from string import rjust
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Open structure file
    mod=loadtxt(home+project_name+'/structure/'+model_name,ndmin=2)
    #Get slip quantities
    iss=2*arange(len(f))
    ids=2*arange(len(f))+1
    ss=sol[iss]
    ds=sol[ids]
    #Get rigidities
    mu=zeros(len(ds))
    trup=zeros(len(ds)*len(rupture_speeds))
    for k in range(len(f)):
        mu[k]=get_mu(mod,f[k,3])
    #Get rupture start times
    for krup in range(len(rupture_speeds)):
        trup[krup*len(ds):(krup+1)*len(ds)]=epi2subfault(epicenter,f,rupture_speeds[krup])
    #Prepare for output
    out=c_[f[:,0:8],ss,ds,f[:,8:10],trup,mu]  #!!!!!!  Have not adjsuted fmt
    outdir=home+project_name+'/output/inverse_models/models/'+run_name+'.'+rjust(str(num),4,'0')+'.inv'
    #CHANGE this to rupture definition as #No  x            y        z(km)      str     dip      rake       rise    dura     slip    ss_len  ds_len rupt_time
    fmtout='%6i\t%.4f\t%.4f\t%6.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\t%.4e%8.1f\t%8.1f\t%8.4f\t%.4e'
    print 'Writing model results to file '+outdir
    savetxt(outdir,out,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
        
    
    
def write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,num):
    '''
    Output synthetics as sac
    '''
    
    from obspy import read
    from numpy import dot,array,savetxt,where,genfromtxt
    from string import rjust
    
    print 'Computing and saving synthetics...'
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
            
        
def write_log(home,project_name,run_name,k,rupture_speeds,next_l,L2,Lm,VR,AIC,Mo,Mw):
    '''
    Write ivnersion sumamry to file
    '''
    from string import rjust
    
    num=rjust(str(k),4,'0')
    f=open(home+project_name+'/output/inverse_models/models/'+run_name+'.'+num+'.log','w')
    f.write('Project: '+project_name+'\n')
    f.write('Run name: '+run_name+'\n')
    f.write('Run number: '+num+'\n')
    f.write('lambda = '+repr(next_l)+'\n')
    f.write('rupture velocities allowed (km/s) = '+str(rupture_speeds)+'\n')
    f.write('L2 = '+repr(L2)+'\n')
    f.write('VR = '+repr(VR)+'\n')
    f.write('Lm = '+repr(Lm)+'\n')
    f.write('AIC = '+repr(AIC)+'\n')
    f.write('M0 = '+repr(Mo)+' N-m\n')
    f.write('Mw = '+repr(Mw)+'\n')
    f.close()
    
    
    
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
    #Check the boundary conditions
    if (top.lower()!='locked' and top.lower()!='free'):
        print 'ERROR: Unknown boundary condition \''+top+'\' at top edge of model'
        return False,False
    if (bottom.lower()!='locked' and bottom.lower()!='free'):
        print 'ERROR: Unknown boundary condition \''+bottom+'\' at bottom edge of model'
        return False,False
    if (right.lower()!='locked' and right.lower()!='free'):
        print 'ERROR: Unknown boundary condition \''+right+'\' at right edge of model'
        return False,False
    if (left.lower()!='locked' and left.lower()!='free'):
        print 'ERROR: Unknown boundary condition \''+left+'\' at left edge of model'
        return False,False
    #No errors, move forward with stencil computation
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
    
    

def get_stats(WG,sol,wd):
    '''
    Compute basic performance metrics of an inversion
    '''
    
    from numpy.linalg import norm
    
    wds=WG.dot(sol)
    L2=norm(wds-wd)
    Lm=norm(sol)
    #Variance reduction
    res=((wd-wds)**2)**0.5
    dnorm=(wd**2)**0.5 #Yes i know this is dumb, shush
    VR=(1-(res.sum()/dnorm.sum()))*100
    AIC=0
    return L2,Lm,VR,AIC
    
    
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
    iss=2*arange(len(f))
    ids=2*arange(len(f))+1
    ss=sol[iss]
    ds=sol[ids]
    #Get total slip
    slip=(ss**2+ds**2)**0.5
    #Get subfault areas in emters
    A=f[:,8]*f[:,9]
    #Get rigidities
    mu=zeros(len(ds))
    for k in range(len(f)):
        mu[k]=get_mu(mod,f[k,3])
    #Compute moments
    M0=mu*A*slip
    #Total up and copute magnitude
    M0=M0.sum()
    Mw=(2./3)*(log10(M0)-9.1)
    
    return M0,Mw