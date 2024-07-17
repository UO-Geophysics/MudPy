'''
D. Melgar 02/2014

Forward modeling routines
'''



def waveforms(home,project_name,rupture_name,station_file,model_name,run_name,integrate,tsunami,hot_start,resample,beta):
    '''
    This routine will take synthetics and apply a slip dsitribution. It will delay each 
    subfault by the appropriate rupture time and linearly superimpose all of them. Output
    will be one sac waveform file per direction of motion (NEU) for each station defined in the
    station_file. Depending on the specified rake angle at each subfault the code will compute 
    the contribution to dip and strike slip directions. It will also compute the moment at that
    subfault and scale it according to the unit amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
       
    OUT:
        Nothing
    '''
    from numpy import loadtxt,genfromtxt,allclose,vstack,deg2rad,array,sin,cos
    from obspy import read,Stream
    import datetime
    import gc
    
    print('Solving for kinematic problem')
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
    stagps_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    #What am I processing v or d ?
    if integrate==1:
        vord='disp'
    else:
        vord='vel'
    #Loop over stations
    for ksta in range(hot_start,len(staname)):
        print('Working on station '+staname[ksta]+' ('+str(ksta+1)+'/'+str(len(staname))+')')
        #Initalize output
        n=Stream()
        e=Stream()
        z=Stream()
        sta=staname[ksta]
        #Loop over sources (Add delays)
        try:
            for k in range(source.shape[0]):
                if k%100==0:
                    print('... working on parameter '+str(k)+' of '+str(len(source)))
                #Get subfault parameters
                nfault='subfault'+str(int(source[k,0])).rjust(4,'0')
                nsub='sub'+str(int(source[k,0])).rjust(4,'0')
                zs=source[k,3]
                ss_slip=source[k,8]
                ds_slip=source[k,9]
                #Rotate
                if beta != None:
                    beta_rot=deg2rad(beta)
                    R=array([[cos(beta_rot),sin(beta_rot)],[-sin(beta_rot),cos(beta_rot)]])
                    rot=R.dot(vstack((ss_slip,ds_slip)))
                    ss_slip=rot[0]
                    ds_slip=rot[1]
                rtime=source[k,12]
                #Where's the data
                strdepth='%.4f' % zs
                if tsunami==False: 
                    syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                else:
                    syn_path=home+project_name+'/GFs/tsunami/'+model_name+'_'+strdepth+'.'+nsub+'/'
                #Get synthetics
                ess=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.e')
                nss=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.n')
                zss=read(syn_path+sta+'.'+nfault+'.SS.'+vord+'.z')
                eds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.e')
                nds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.n')
                zds=read(syn_path+sta+'.'+nfault+'.DS.'+vord+'.z')
                #Decide if resampling is required
                if resample!=None:
                    if resample < (1/ess[0].stats.delta): #Downsample
                        ess[0].resample(resample)
                        nss[0].resample(resample)
                        zss[0].resample(resample)
                        eds[0].resample(resample)
                        nds[0].resample(resample)
                        zds[0].resample(resample)
                    elif resample > (1/ess[0].stats.delta): #Upsample
                        upsample(ess,1./resample)
                        upsample(nss,1./resample)
                        upsample(zss,1./resample)
                        upsample(eds,1./resample)
                        upsample(nds,1./resample)
                        upsample(zds,1./resample)
                dt=ess[0].stats.delta
                #Time shift them according to subfault rupture time
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
                if allclose((ss_slip**2+ds_slip**2)**0.5,0)==False:  #Only add things that matter
                    log=log+nfault+', SS='+str(ss_slip)+', DS='+str(ds_slip)+'\n'
                    #A'ight, add 'em up
                    etotal=add_traces(ess,eds,ss_slip,ds_slip)
                    ntotal=add_traces(nss,nds,ss_slip,ds_slip)
                    ztotal=add_traces(zss,zds,ss_slip,ds_slip)
                    #Add to previous subfault's results
                    e=add_traces(e,etotal,1,1)
                    n=add_traces(n,ntotal,1,1)
                    z=add_traces(z,ztotal,1,1)
                else:
                    log=log+"No slip on subfault "+nfault+', ignoring it...\n'
                gc.collect()
            #Save results
            e.write(outpath+run_name+'.'+sta+'.'+vord+'.e',format='SAC')
            n.write(outpath+run_name+'.'+sta+'.'+vord+'.n',format='SAC')
            z.write(outpath+run_name+'.'+sta+'.'+vord+'.u',format='SAC')
        except:
            print('An error coccured, skipping station')
    f=open(logpath+'waveforms.'+now+'.log','a')
    f.write(log)
    f.close()
        
        
        
def waveforms_matrix(home,project_name,fault_name,rupture_name,station_file,GF_list,
                model_name,run_name,epicenter,time_epi,integrate,tsunami,hot_start,
                resample,beta,rupture_speed,num_windows,dt,NFFT):
    '''
    To supplant waveforms() it needs to include resmapling and all that jazz...
    
    This routine will take synthetics and apply a slip dsitribution. It will delay each 
    subfault by the appropriate rupture time and linearly superimpose all of them. Output
    will be one sac waveform file per direction of motion (NEU) for each station defined in the
    station_file. Depending on the specified rake angle at each subfault the code will compute 
    the contribution to dip and strike slip directions. It will also compute the moment at that
    subfault and scale it according to the unit amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
       
    OUT:
        Nothing
    '''
    from numpy import loadtxt,genfromtxt,allclose,vstack,deg2rad,array,sin,cos,where,zeros,arange
    from obspy import read,Stream,Trace
    import datetime
    import gc
    from mudpy.inverse import getG
    from linecache import getline
    from os import remove
    
    print('Solving for kinematic problem')
    #Output where?
    outpath=home+project_name+'/output/forward_models/'
    logpath=home+project_name+'/logs/'
    log=''
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #load source
    #source=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    #Load stations
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    #Now load rupture model
    mss=genfromtxt(home+project_name+'/forward_models/'+rupture_name,usecols=8)
    mds=genfromtxt(home+project_name+'/forward_models/'+rupture_name,usecols=9)
    m=zeros(2*len(mss))
    i=arange(0,2*len(mss),2)
    m[i]=mss
    i=arange(1,2*len(mds),2)
    m[i]=mds
    #What am I processing v or d ?
    if integrate==1:
        vord='disp'
    else:
        vord='vel'
    #Load gflist
    gfsta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,skip_header=1,dtype='U')
    #Loop over stations
    for ksta in range(hot_start,len(staname)):
        print('Working on station '+staname[ksta]+' ('+str(ksta+1)+'/'+str(len(staname))+')')
        #Initalize output
        n=Stream()
        e=Stream()
        z=Stream()
        sta=staname[ksta]
        #Make dummy data
        ndummy=Stream(Trace())
        ndummy[0].data=zeros(int(NFFT))
        ndummy[0].stats.delta=dt
        ndummy[0].stats.starttime=time_epi
        edummy=ndummy.copy()
        udummy=ndummy.copy()
        ndummy.write(home+project_name+'/data/waveforms/'+sta+'.'+vord+'.n',format='SAC')
        edummy.write(home+project_name+'/data/waveforms/'+sta+'.'+vord+'.e',format='SAC')
        udummy.write(home+project_name+'/data/waveforms/'+sta+'.'+vord+'.u',format='SAC')
        #Extract only one station from GF_list file
        ista=int(where(gfsta==staname[ksta])[0])+2
        #Make mini GF_file
        tmpgf='tmpfwd.gflist'
        try:
            remove(home+project_name+'/data/station_info/'+tmpgf)
        except:
            pass
        gflist=getline(home+project_name+'/data/station_info/'+GF_list,ista)
        f=open(home+project_name+'/data/station_info/'+tmpgf,'w')
        f.write('# Headers\n')
        f.write(gflist)
        f.close()
        #Get matrix for one station to all sources
        G_from_file=False
        G_name='tmpfwd'
        G=getG(home,project_name,fault_name,model_name,tmpgf,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,decimate=None,bandpass=None,tsunami=False)
        # Matrix multiply and separate data streams
        d=G.dot(m)
        n=ndummy.copy()
        e=edummy.copy()
        u=udummy.copy()
        ncut=len(d)/3
        n[0].data=d[0:ncut]
        e[0].data=d[ncut:2*ncut]
        u[0].data=d[2*ncut:3*ncut]
        # Write to file
        n.write(home+project_name+'/output/forward_models/'+run_name+'.'+gfsta[ksta]+'.'+vord+'.n',format='SAC')
        e.write(home+project_name+'/output/forward_models/'+run_name+'.'+gfsta[ksta]+'.'+vord+'.e',format='SAC')
        u.write(home+project_name+'/output/forward_models/'+run_name+'.'+gfsta[ksta]+'.'+vord+'.u',format='SAC')


def waveforms_fakequakes(home,project_name,fault_name,rupture_list,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function='dreger',zeta=0.2,
                stf_falloff_rate=4.0,rupture_name=None,epicenter=None,time_epi=None,hot_start=0,
                ncpus=1):
    '''
    To supplant waveforms_matrix() it needs to include resmapling and all that jazz...
    
    Instead of doing matrix multiplication one stations at a time, do it for all stations
    
    This routine will take synthetics and apply a slip dsitribution. It will delay each 
    subfault by the appropriate rupture time and linearly superimpose all of them. Output
    will be one sac waveform file per direction of motion (NEU) for each station defined in the
    station_file. Depending on the specified rake angle at each subfault the code will compute 
    the contribution to dip and strike slip directions. It will also compute the moment at that
    subfault and scale it according to the unit amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
       
    OUT:
        Nothing
    '''
    from numpy import genfromtxt,array
    from obspy import UTCDateTime
    import datetime
    import time
    from dask.distributed import Client
    
    print('Solving for kinematic problem(s)')
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    
    #load source names
    if rupture_list==None:
        #all_sources=array([home+project_name+'/forward_models/'+rupture_name])   
        all_sources=array([rupture_name])  
    else:
        all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
        all_sources = array(all_sources, ndmin=1)  # in case only 1 entry

    
    
    
    #Load all synthetics
    print('... loading all synthetics into memory')
    Nss,Ess,Zss,Nds,Eds,Zds=load_fakequakes_synthetics(home,project_name,fault_name,model_name,GF_list,
                                                       G_from_file,G_name)
    print('... ... ... done')
    

    # Need to know how many sites and delta t
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    Nsta=len(staname)
    dt = Nss[0].stats.delta
    
    # Need epicetrnal time from log file to trim synthetics
    print('... reading epicentral time from log file. REMEMBER: All ruptures in ruptures.list should have a common epicentral time')
    log_file = home + project_name + '/output/ruptures/' + all_sources[0].replace('.rupt','.log')
    flog = open(log_file,'r')
    while True:
        line = flog.readline()
        if 'Hypocenter time' in line:                
            hypocenter_time = line.replace('Hypocenter time: ','')
            hypocenter_time = UTCDateTime(hypocenter_time)
            break
        elif line == '':
            break 
    flog.close()
    
    #Now get the impulse response G for all sites and all subfaults
    print('... broadcast to G matrix')
    Gimpulse_all = stream2matrix(Nss,Ess,Zss,Nds,Eds,Zds,hypocenter_time,Nsta)
    print('... ... done')
    
    #clean up
    print('... clean up: removing obspy streams from memory')
    del Nss,Ess,Zss,Nds,Eds,Zds
    print('... ... done')
    
    # print('... opening DASK client')
    # client = Client(n_workers=ncpus)
    
    
    #Now loop over rupture models
    for ksource in range(hot_start,len(all_sources)):
       
        rupture_name = all_sources[ksource]
        print('... solving for source '+str(ksource)+' of '+str(len(all_sources))+': '+rupture_name)
        
        if rupture_list != None:
            #Get epicentral time
            epicenter,time_epi = read_fakequakes_hypo_time(home,project_name,rupture_name)
            forward = False
        else:
            forward = True #This controls where we look for the rupture file
        
        # Put in matrix
        t1=time.time()
        # m,G=get_fakequakes_G_and_m(Nss,Ess,Zss,Nds,Eds,Zds,home,project_name,rupture_name,time_epi,GF_list,epicenter,NFFT,source_time_function,stf_falloff_rate,zeta,forward=forward)
        m,G = get_fakequakes_G_and_m(Gimpulse_all,home,project_name,rupture_name,time_epi,GF_list,
                                   epicenter,NFFT,source_time_function,stf_falloff_rate,zeta=zeta,
                                   forward=forward,dt=dt,ncpus=ncpus,reconvolution=False)
        t2=time.time()
        print('... ... slip rate convolutions completed: wall time {:.1f}s'.format(t2-t1))
        # Solve
        print('... ... solve for waveforms (ncpus = %d)' % ncpus)
        waveforms = G.dot(m)
        #Write output
        write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,GF_list,NFFT,time_epi,dt)
        print('... ... finished with this rupture')
        
    # print('... clean up: closing distributed client')
    # client.close()
        
   
    
   
def stream2matrix(Nss,Ess,Zss,Nds,Eds,Zds,hypocenter_time,Nsta):
    '''
    
    Converts stream objects to a properly formatted G matrix

    Parameters
    ----------
    Nss,Ess ... : Stream objects with synthetics for each component of motion and
                    for the ss and ds rake angles


    hypocenter_time : UTCDateTime object with hypocentral timne, it's used to troim
        the syntehtics to a single, common time

    Nsta : int, number of stations being processed


    Returns
    -------
    G = [ Nss_sta1_sf1 Nds sta1_ sf1     Nss_sta1_sf2 Nds sta1_ sf2 ...
          Ess_sta1_sf1 Eds sta1_ sf1     Nss_sta1_sf3 Nds sta1_ sf3 ...
          Zss_sta1_sf1 Zds sta1_ sf1     Nss_sta1_sf3 Nds sta1_ sf3 ...
          Nss_sta2_sf1 Nds sta2_ sf1
          Ess_sta2_sf1 Eds sta2_ sf1
          Zss_sta2_sf1 Zds sta2_ sf1
          ...

    '''    
    
    from numpy import ones
    
    #How many time points in snthethics
    Npts = Nss[0].stats.npts
    
    #How many sources?
    Nsources = int(len(Nss) / Nsta)
    
    #Initalize G matrix
    G = ones((Npts*Nsta*3,2*Nsources))
    
    k=0
    row_start = 0
    row_end = Npts
    for ksta in range(Nsta):
        for ksub in range(Nsources):
            

            
            #Assign individual GFs
            nss =  Nss[k]
            ess =  Ess[k]
            zss =  Zss[k]
            
            nds =  Nds[k]
            eds =  Eds[k]
            zds =  Zds[k]
            
            # Trim to hypocentral time and pad
            nss = trim_to_hypo_time_and_pad(nss,hypocenter_time,Npts)
            ess = trim_to_hypo_time_and_pad(ess,hypocenter_time,Npts)
            zss = trim_to_hypo_time_and_pad(zss,hypocenter_time,Npts)

            nds = trim_to_hypo_time_and_pad(nds,hypocenter_time,Npts)
            eds = trim_to_hypo_time_and_pad(eds,hypocenter_time,Npts)
            zds = trim_to_hypo_time_and_pad(zds,hypocenter_time,Npts)
            
            # Done, now place in correct position in matrix
            G[row_start:row_end,2*ksub] = nss
            G[row_start:row_end,2*ksub+1] = nds
            
            G[row_start+Npts:row_end+Npts,2*ksub] = ess
            G[row_start+Npts:row_end+Npts,2*ksub+1] = eds
    
            G[row_start+2*Npts:row_end+2*Npts,2*ksub] = zss
            G[row_start+2*Npts:row_end+2*Npts,2*ksub+1] = zds
            
            k += 1
    
        row_start += 3*Npts
        row_end += 3*Npts
    
    return G
        
        

  
        
        
def tele_waveforms(home,project_name,fault_name,rupture_list,GF_list_teleseismic,
                model_name,run_name,G_from_file,G_name,source_time_function='dreger',zeta=0.2,
                stf_falloff_rate=4.0,rupture_name=None,epicenter=None,time_epi=None,
                hot_start=0,decimation_factor=1,reconvolution=True,ncpus=1):
    '''

    '''
    from numpy import genfromtxt,array
    import datetime
    import gc
    from obspy import UTCDateTime
    
    print('Solving for kinematic problem(s)')
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    
    #load source names
    if rupture_list==None: 
        all_sources=array([rupture_name])  
    else:
        all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
    
    #Load all synthetics
    print('... loading all synthetics into memory')
    Nss,Ess,Zss,Nds,Eds,Zds=load_fakequakes_tele_synthetics(home,project_name,fault_name,model_name,GF_list_teleseismic,G_from_file,G_name,decimation_factor)
    print('... ... done')
    
    
    # Need to know how many sites and delta t
    station_file=home+project_name+'/data/station_info/'+GF_list_teleseismic
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    Nsta=len(staname)
    Npoints=Nss[0].stats.npts
    dt=Nss[0].stats.delta
    
    # Need epicetrnal time from log file to trim synthetics
    print('... reading epicentral time from log file. REMEMBER: All ruptures in ruptures.list should have a common epicentral time')
    log_file = home + project_name + '/output/ruptures/' + all_sources[0].replace('.rupt','.log')
    flog = open(log_file,'r')
    while True:
        line = flog.readline()
        if 'Hypocenter time' in line:                
            hypocenter_time = line.replace('Hypocenter time: ','')
            hypocenter_time = UTCDateTime(hypocenter_time)
            break
        elif line == '':
            break 
    flog.close()
    
    #Now get the impulse response G for all sites and all subfaults
    print('... broadcast to G matrix')
    Gimpulse_all = stream2matrix(Nss,Ess,Zss,Nds,Eds,Zds,hypocenter_time,Nsta)
    print('... ... done')
    
    #clean up
    print('... clean up: removing obspy streams from memory')
    del Nss,Ess,Zss,Nds,Eds,Zds
    print('... ... done')
    
    
    #Now loop over rupture models
    for ksource in range(hot_start,len(all_sources)):
        print('... solving for source '+str(ksource)+' of '+str(len(all_sources)))
        rupture_name=all_sources[ksource]
        print('... ... '+rupture_name)
        
        if rupture_list!=None:
            #Get epicentral time
            epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
            forward=False
        else:
            forward=True #This controls where we look for the rupture file
        
        # Put in matrix
        if source_time_function == 'gauss_prem_i2s':
            reconvolution = True
        else:
            reconvolution = False
        
        m,G = get_fakequakes_G_and_m(Gimpulse_all,home,project_name,rupture_name,time_epi,GF_list_teleseismic,
                                   epicenter,Npoints,source_time_function,stf_falloff_rate,zeta=zeta,
                                   forward=forward,dt=dt,ncpus=ncpus,reconvolution=reconvolution)
        # Solve
        waveforms=G.dot(m)
        print('... ... DONE! saving waveforms')
        #Write output
        write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,GF_list_teleseismic,Npoints,time_epi,dt)
        
        #Delete variables and garbage collect
        del G
        gc.collect()
        


def waveforms_fakequakes_dynGF(home,project_name,fault_name,rupture_list,GF_list,dynamic_GFlist,dist_threshold,
                         model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function='dreger',
                         stf_falloff_rate=4.0,rupture_name=None,epicenter=None,time_epi=None,
                         hot_start=0,ncpus=1):
    '''
    To supplant waveforms_matrix() it needs to include resmapling and all that jazz...
    
    Instead of doing matrix multiplication one stations at a time, do it for all stations
    
    This routine will take synthetics and apply a slip dsitribution. It will delay each 
    subfault by the appropriate rupture time and linearly superimpose all of them. Output
    will be one sac waveform file per direction of motion (NEU) for each station defined in the
    station_file. Depending on the specified rake angle at each subfault the code will compute 
    the contribution to dip and strike slip directions. It will also compute the moment at that
    subfault and scale it according to the unit amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to de displacement
        dynamic_GFlist: A boolean (Useless for current version)
        dist_threshold: A float, station to the closest subfault must be closer to this distance, otherwise skip it.
        
       
    OUT:
        Nothing
    '''
    from numpy import genfromtxt,array
    import datetime
    import numpy as np
    try:
        import multiprocessing as mp
        import time
        use_parallel=True
    except:
        print('Parallel waveform generation is unavailable...')
        print('Please pip install multiprocessing')
        print('Otherwise, this is now continue go with single cpu')
        use_parallel=False
    
    print('Solving for kinematic problem(s)')
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')

    #load source names
    if rupture_list==None:
        #all_sources=array([home+project_name+'/forward_models/'+rupture_name])   
        all_sources=array([rupture_name])
    else:
        all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')

    #Load all synthetics
    print('... loading all synthetics into memory')
    Nss,Ess,Zss,Nds,Eds,Zds=load_fakequakes_synthetics(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name)
    print('... ... done loading synthetics')

    #Load GF_list
    STAname=np.genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[0],skip_header=1,dtype='S6')
    STA={sta.decode():nsta for nsta,sta in enumerate(STAname)}

    def loop_sources(ksource):
        print('...Solving for source '+str(ksource)+' of '+str(len(all_sources)))
        rupture_name=all_sources[ksource]
        ###make new GF_list(dynamic GF_list)###
        new_GF_list='subGF_list.rupt%06d.gflist'%(ksource) #make new GF_list for only close stations
        sta_close_to_rupt(home,project_name,rupture_name,GF_list,dist_threshold,new_GF_list)
        
        if rupture_list!=None:
            #Get epicentral time
            epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
            forward=False
        else:
            forward=True #This controls where we look for the rupture file
        
        # Put in matrix
        #use only close stations (dynamic GF_list)
        m,G=get_fakequakes_G_and_m_dynGF(Nss,Ess,Zss,Nds,Eds,Zds,home,project_name,rupture_name,time_epi,STA,new_GF_list,epicenter,NFFT,source_time_function,stf_falloff_rate,forward=forward)

        # Solve
        waveforms=G.dot(m)
        ##Write output
        write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,new_GF_list,NFFT,time_epi,dt)

    if use_parallel:
        print('Using cpu=',ncpus)
        #Parallel(n_jobs=ncpus,backend='loky')(delayed(loop_sources)(ksource) for ksource in range(hot_start,len(all_sources))) #Oh no, this is REALLY slow why?
        ps=[]
        for ksource in range(hot_start,len(all_sources)):
            p = mp.Process(target=loop_sources, args=([ksource]) )
            ps.append(p)
        #start running
        queue=np.zeros(len(ps))
        for i_ps in range(len(ps)):
            if (i_ps%10==0):
                print('now at',i_ps,'out of',len(ps))
            ps[i_ps].start()
            queue[i_ps]=1
            while True:
                #check running status
                running_idx=np.where(queue==1)[0]
                for ri in running_idx:
                    if not(ps[ri].is_alive()):
                        #finnish running, close
                        ps[ri].join()
                        #ps[ri]=0 #.join() still has issue...
                        queue[ri]=0
                    else:
                        continue
                if len(np.where(queue==1)[0])<=ncpus:
                    #add a process
                    #print('number of processer=',len(np.where(queue==1)[0]),'add a process,now at',i)
                    break
                else:
                    #print('number of queue reaches max:',nprocess,'try again later,now at',i)
                    time.sleep(0.5) #wait and try again later
        #Final check if all the processes are done
        while True:
            if np.sum([nps.is_alive() for nps in ps ])==0:
                break
            else:
                time.sleep(1)

    else:
        #Now loop over rupture models
        for ksource in range(hot_start,len(all_sources)):
            loop_sources(ksource)
            '''
            print('... solving for source '+str(ksource)+' of '+str(len(all_sources)))
            rupture_name=all_sources[ksource]
            print(rupture_name)
            
            ###make new GF_list(dynamic GF_list)###
            new_GF_list='subGF_list.rupt%06d.gflist'%(ksource) #make new GF_list for only close stations
            sta_close_to_rupt(home,project_name,rupture_name,GF_list,dist_threshold,new_GF_list)
            
            if rupture_list!=None:
                #Get epicentral time
                epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
                forward=False
            else:
                forward=True #This controls where we look for the rupture file
                
            # Put in matrix
            #use only close stations (dynamic GF_list)
            m,G=get_fakequakes_G_and_m_dynGF(Nss,Ess,Zss,Nds,Eds,Zds,home,project_name,rupture_name,time_epi,STA,new_GF_list,epicenter,NFFT,source_time_function,stf_falloff_rate,forward=forward)
            
            # Solve
            waveforms=G.dot(m)
            ##Write output
            #write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,GF_list,NFFT,time_epi,dt)
            #Write output
            write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,new_GF_list,NFFT,time_epi,dt)
            '''





def hf_waveforms(home,project_name,fault_name,rupture_list,GF_list,model_name,run_name,dt,NFFT,G_from_file,
            G_name,rise_time_depths,moho_depth_in_km,ncpus,source_time_function='dreger',duration=100.0,
            stf_falloff_rate=4.0,hf_dt=0.02,Pwave=False,Swave=True,hot_start=0,stress_parameter=50,
            high_stress_depth=1e4,kappa=None,Qexp=0.6,Qmethod='shallowest',scattering='off',Qc_exp=0,
            baseline_Qc=100):

    '''
    Make semistochastic high frequency accelerograms
    '''
    from numpy import genfromtxt,ones
    from mudpy import hfsims
    
    
    all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
    #Now loop over rupture models
    Nsources=all_sources.size
    for ksource in range(hot_start,Nsources):        
        print('... solving HF waveforms for source '+str(ksource)+' of '+str(Nsources))
        if Nsources>1:
            rupture_name=all_sources[ksource]
        else:
            rupture_name=str(all_sources)
        
        #Get station info from GF_list
        sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[0],dtype='U')
        lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
        sta_lon=lonlat[:,0]
        sta_lat=lonlat[:,1]
        
        #Multi kappa?
        if kappa == 'gflist': #from station files
            station_kappa = genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[-1])
        elif kappa == None: #default values
            station_kappa = 0.04*ones(len(sta))
        else: #one single value for all sites that is NOT the default
            station_kappa = kappa*ones(len(sta))
        
        comp=['N','E','Z']  #components to loop through
        
        #Now loop over stations
        for ksta in range(len(sta)):
            #Now loop over components N,E,Z
            
            kappa = station_kappa[ksta] # assign station specific kappa values
            
            for kcomp in range(len(comp)):
                #HF_sims stochastic simulation for single station, component
                make_parallel_hfsims(home,project_name,rupture_name,ncpus,sta[ksta],sta_lon[ksta],sta_lat[ksta],
                    comp[kcomp],model_name,rise_time_depths[0],rise_time_depths[1],moho_depth_in_km,total_duration=duration,hf_dt=hf_dt,
                    Pwave=Pwave,Swave=Swave,stress_parameter=stress_parameter,high_stress_depth=high_stress_depth,Qexp=Qexp,
                    Qmethod=Qmethod,scattering=scattering,Qc_exp=Qc_exp,baseline_Qc=baseline_Qc)

                #Combine the separate MPI outputs into one full waveform
                write_parallel_hfsims(home,project_name,rupture_name,sta[ksta],comp[kcomp],remove=True)



def make_parallel_hfsims(home,project_name,rupture_name,ncpus,sta,sta_lon,sta_lat,component,model_name,rise_time_depths0,
                         rise_time_depths1,moho_depth_in_km,total_duration,hf_dt,Pwave,Swave,stress_parameter,
                         kappa=0.04,Qexp=0.6,high_stress_depth=30,Qmethod='shallowest',scattering='off',Qc_exp=0,
                         baseline_Qc=100):


    '''
    Set up for MPI calculation of HF stochastics
    '''
    from numpy import savetxt,arange,genfromtxt,where
    from os import environ
    import subprocess
    from shlex import split
    
    #Calculate the necessary full-fault parameters before splitting up the faults over your ncpus
    rupture=home+project_name+'/output/ruptures/'+rupture_name
    fault=genfromtxt(rupture)
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    subfault_M0=slip*fault[:,10]*fault[:,11]*fault[:,13]
    subfault_M0=subfault_M0*1e7 #to dyne-cm
    M0=subfault_M0.sum()
    i=where(slip>0)[0]
    N=len(i)
    #Split rupture file into ncpu rupture files
    for k in range(ncpus):
        i=arange(k,len(fault),ncpus)
        mpi_source=fault[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6E'
        savetxt(home+project_name+'/output/ruptures/mpi_rupt.'+str(k)+'.'+rupture_name,mpi_source,fmt=fmt)
    #Make mpi system call
    print("MPI: Starting Stochastic High Frequency Simulation on ", ncpus, "CPUs")
    mud_source=environ['MUD']+'/src/python/mudpy/'
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'hfsims_parallel.py run_parallel_hfsims '+home+' '+project_name+' '+rupture_name+' '+str(N)+' '+str(M0)+' '+sta+' '+str(sta_lon)+' '+str(sta_lat)+' '+model_name+' '+str(rise_time_depths0)+' '+str(rise_time_depths1)+' '+str(moho_depth_in_km)+' '+component+' '+str(total_duration)+' '+str(hf_dt)+' '+str(stress_parameter)+' '+str(kappa)+' '+str(Qexp)+' '+str(Pwave)+' '+str(Swave)+' '+str(high_stress_depth)+' '+str(Qmethod)+' '+str(scattering)+' '+str(Qc_exp)+' '+str(baseline_Qc)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
    
    
                
def run_hf_waveforms(home,project_name,fault_name,rupture_list,GF_list,model_name,run_name,dt,NFFT,G_from_file,
            G_name,rise_time_depths,moho_depth_in_km,source_time_function='dreger',duration=100.0,
            stf_falloff_rate=4.0,hf_dt=0.02,Pwave=False,hot_start=0,stress_parameter=50,
            high_stress_depth=1e4,Qexp=0.6):
    '''
    Make semistochastic high frequency accelerograms
    '''                        
    from numpy import genfromtxt
    #import datetime
    from mudpy import hfsims
    
    #load list of source names
    all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')

    #Now loop over rupture models
    Nsources=all_sources.size
    for ksource in range(hot_start,Nsources):
        
        print('... solving HF waveforms for source '+str(ksource)+' of '+str(Nsources))
        if Nsources>1:
            rupture_name=all_sources[ksource]
        else:
            rupture_name=str(all_sources)
        print(rupture_name)
        
        #Get epicentral time
        epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
        
        #Get station info from GF_list
        sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[0],dtype='U')
        lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
        
        comp=['N','E','Z']  #components to loop through
        
        #Now loop over stations
        for ksta in range(len(sta)):
            #Now loop over components N,E,Z
            for kcomp in range(len(comp)):
                #HF_sims stochastic simulation for single station, component
                hf_waveform=hfsims.stochastic_simulation(home,project_name,rupture_name,sta[ksta],lonlat[ksta,:],time_epi,
                        model_name,rise_time_depths,moho_depth_in_km,hf_dt=hf_dt,total_duration=duration,component=comp[kcomp],
                        Pwave=Pwave,stress_parameter=stress_parameter,high_stress_depth=high_stress_depth)
                #Write to file
                write_fakequakes_hf_waveforms_one_by_one(home,project_name,rupture_name,hf_waveform,comp[kcomp])
    


def write_parallel_hfsims(home,project_name,rupture_name,station,component,remove=False):
    '''
    Combine outputs from parallel hfsims into one complete waveform.
    Delete the parallel hfsims outputs when complete.
    '''
    from glob import glob
    from obspy import read
    import os
    
    #rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new
    rupture = rupture_name.rsplit('.', 1)[0]
    parallel_waveforms=glob(home+project_name+'/output/waveforms/'+rupture+'/'+station+'.HN'+component+'.[0-9][0-9][0-9].sac')
    #print "Number of MPI outputs: " + str(len(parallel_waveforms))
    if len(parallel_waveforms)==0:
        print("No waveforms at this location to add in")
    else:
        for r in range(len(parallel_waveforms)):
            #print parallel_waveforms[r]
            st=read(parallel_waveforms[r])
            if r==0:
                #stdout.write('      [.'); stdout.flush()
                complete_waveform=st.copy()
            else:
                #stdout.write('.'); stdout.flush()
                complete_waveform[0].data=complete_waveform[0].data+st[0].data
        complete_waveform.write(home+project_name+'/output/waveforms/'+rupture+'/'+station+'.HN'+component+'.mpi.sac',format='SAC')
        if remove==True:
            for r in range(len(parallel_waveforms)):
                os.remove(parallel_waveforms[r])
        else:
            print("Keeping all parallel output files because you didn't tell me to delete them")
                       
def write_fakequakes_hf_waveforms_one_by_one(home,project_name,rupture_name,hf_trace,component):
    '''
    write HF waveforms to file as they are created, station by station, component by component
    '''

    from os import path,makedirs
    
    #Where am I writing to?
    #rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new
    rupture = rupture_name.rsplit('.', 1)[0]
    directory=home+project_name+'/output/waveforms/'+rupture+'/'
    
    #Check if dir exists if not then create it
    if not path.exists(directory):
        makedirs(directory)   
        
    sta=hf_trace.stats.station
    hf_trace.write(directory+sta+'.HN'+component+'.sac',format='SAC')


def write_fakequakes_hf_waveforms(home,project_name,rupture_name,n,e,z):
    '''
    write HF waveforms to file
    '''

    from os import path,makedirs
    from numpy import genfromtxt,squeeze,sqrt
    from obspy import Stream,Trace
    
    #Where am I writing to?
    #rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new
    rupture = rupture_name.rsplit('.', 1)[0]
    directory=home+project_name+'/output/waveforms/'+rupture+'/'
    
    #Check if dir exists if not then create it
    if not path.exists(directory):
        makedirs(directory)   
        
    for ksta in range(len(n)):
        sta=n[ksta].stats.station
        n[ksta].write(directory+sta+'.HNN.sac',format='SAC')
        e[ksta].write(directory+sta+'.HNE.sac',format='SAC')
        z[ksta].write(directory+sta+'.HNZ.sac',format='SAC')
        
                                                                                      

def match_filter(home,project_name,fault_name,rupture_list,GF_list,
        zero_phase=False,order=2,fcorner_low=1.0, fcorner_high=1.0):
    '''
    match filter waveforms
    '''
    
    from numpy import genfromtxt,where,r_,diff,interp
    from obspy import read
    
    #read stations list
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U')
    ista=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=4)
    ista=where(ista==1)[0]
    sta=sta[ista]
    
    #read ruptures
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
    Nruptures=ruptures.size
    
    for krup in range(Nruptures):
        
        if Nruptures>1:
            rupture_name=ruptures[krup]
        else:
            rupture_name=str(ruptures)
        rupture=rupture_name.replace('.rupt','')
        directory=home+project_name+'/output/waveforms/'+rupture+'/'
        
        print('Running matched filter for all stations for rupture '+ rupture_name)
        
        #Get epicentral time
        epicenter,time_epi=read_fakequakes_hypo_time(home,project_name,rupture_name)
        
        for ksta in range(len(sta)):
            
            #Read HF accelerograms
            try:
                hf_n=read(directory+sta[ksta]+'.HNN.sac')
                hf_e=read(directory+sta[ksta]+'.HNE.sac')
                hf_z=read(directory+sta[ksta]+'.HNZ.sac')
            except:
                hf_n=read(directory+sta[ksta]+'.HNN.mpi.sac')
                hf_e=read(directory+sta[ksta]+'.HNE.mpi.sac')
                hf_z=read(directory+sta[ksta]+'.HNZ.mpi.sac')
        
            #Read LF displacements
            lf_n=read(directory+sta[ksta]+'.LYN.sac')
            lf_e=read(directory+sta[ksta]+'.LYE.sac')
            lf_z=read(directory+sta[ksta]+'.LYZ.sac')
            
            
            #Diff LF to acceleration
            dt=lf_n[0].stats.delta
            
            lf_n[0].data=r_[0,diff(lf_n[0].data)/dt]
            lf_e[0].data=r_[0,diff(lf_e[0].data)/dt]
            lf_z[0].data=r_[0,diff(lf_z[0].data)/dt]
            
            lf_n[0].data=r_[0,diff(lf_n[0].data)/dt]
            lf_e[0].data=r_[0,diff(lf_e[0].data)/dt]
            lf_z[0].data=r_[0,diff(lf_z[0].data)/dt]
            
            #Apply filters
            fsample=1./hf_n[0].stats.delta
            hf_n[0].data=highpass(hf_n[0].data,fcorner_high,fsample,order,zerophase=zero_phase)
            hf_e[0].data=highpass(hf_e[0].data,fcorner_high,fsample,order,zerophase=zero_phase)
            hf_z[0].data=highpass(hf_z[0].data,fcorner_high,fsample,order,zerophase=zero_phase)

            fsample=1./lf_n[0].stats.delta
            lf_n[0].data=lowpass(lf_n[0].data,fcorner_low,fsample,order,zerophase=zero_phase)
            lf_e[0].data=lowpass(lf_e[0].data,fcorner_low,fsample,order,zerophase=zero_phase)
            lf_z[0].data=lowpass(lf_z[0].data,fcorner_low,fsample,order,zerophase=zero_phase)
            
            #Resample LF to HF sample rate
            tinterp=hf_n[0].times() #interpolation time vector
            lf_n[0].data=interp(tinterp,lf_n[0].times(),lf_n[0].data)
            lf_e[0].data=interp(tinterp,lf_e[0].times(),lf_e[0].data)
            lf_z[0].data=interp(tinterp,lf_z[0].times(),lf_z[0].data)
            lf_n[0].stats.delta=hf_n[0].stats.delta
            lf_e[0].stats.delta=hf_n[0].stats.delta
            lf_z[0].stats.delta=hf_n[0].stats.delta
            
            #Trim HF and LF to match each otehrs start and end times
            
            if lf_n[0].stats.endtime < hf_n[0].stats.endtime:
                tend=lf_n[0].stats.endtime
            else:
                tend=hf_n[0].stats.endtime
            
            lf_n[0].trim(starttime=time_epi,endtime=tend)
            hf_n[0].trim(starttime=time_epi,endtime=tend)
            lf_e[0].trim(starttime=time_epi,endtime=tend)
            hf_e[0].trim(starttime=time_epi,endtime=tend)
            lf_z[0].trim(starttime=time_epi,endtime=tend)
            hf_z[0].trim(starttime=time_epi,endtime=tend)
            
            #Add together for final waveform
            bb_n=hf_n.copy()
            bb_n[0].data=lf_n[0].data+hf_n[0].data
            bb_e=hf_e.copy()
            bb_e[0].data=lf_e[0].data+hf_e[0].data
            bb_z=hf_z.copy()
            bb_z[0].data=lf_z[0].data+hf_z[0].data
            
            #Write to file
            bb_n.write(directory+sta[ksta]+'.bb.HNN.sac',format='SAC')
            bb_e.write(directory+sta[ksta]+'.bb.HNE.sac',format='SAC')
            bb_z.write(directory+sta[ksta]+'.bb.HNZ.sac',format='SAC')
            
            
            
            

def load_fakequakes_synthetics(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name):
    ''''
    Load the miniseed files with all the synthetics
    '''
    from numpy import genfromtxt,loadtxt
    from obspy import read    

    vord='disp'
    if G_from_file==True: #load from file
        print('... ... read from miniSEED into stream objects')
        Eds=read(home+project_name+'/GFs/matrices/'+G_name+'.Eds.'+vord+'.mseed')
        Nds=read(home+project_name+'/GFs/matrices/'+G_name+'.Nds.'+vord+'.mseed')
        Zds=read(home+project_name+'/GFs/matrices/'+G_name+'.Zds.'+vord+'.mseed')
        Ess=read(home+project_name+'/GFs/matrices/'+G_name+'.Ess.'+vord+'.mseed')
        Nss=read(home+project_name+'/GFs/matrices/'+G_name+'.Nss.'+vord+'.mseed')
        Zss=read(home+project_name+'/GFs/matrices/'+G_name+'.Zss.'+vord+'.mseed')
    else: #assemble G one data type at a time, just displacememnt right now
        #Load station info
        station_file=home+project_name+'/data/station_info/'+GF_list
        staname=genfromtxt(station_file,dtype="U",usecols=0)
        Nsta=len(staname)
        #Load fault model
        source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
        Nfaults=source.shape[0] #Number of subfaults
        kindex=0
        for ksta in range(Nsta):
            print('... ... reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
            for kfault in range(Nfaults):
                #Get subfault GF directory
                nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                strdepth='%.4f' % source[kfault,3]
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
                kindex+=1
        print('... done, writting synthetics to miniSEED, hang on this might take a minute or two.')
        Ess.write(home+project_name+'/GFs/matrices/'+G_name+'.Ess.'+vord+'.mseed',format='MSEED')
        Nss.write(home+project_name+'/GFs/matrices/'+G_name+'.Nss.'+vord+'.mseed',format='MSEED')
        Zss.write(home+project_name+'/GFs/matrices/'+G_name+'.Zss.'+vord+'.mseed',format='MSEED')
        Eds.write(home+project_name+'/GFs/matrices/'+G_name+'.Eds.'+vord+'.mseed',format='MSEED')
        Nds.write(home+project_name+'/GFs/matrices/'+G_name+'.Nds.'+vord+'.mseed',format='MSEED')
        Zds.write(home+project_name+'/GFs/matrices/'+G_name+'.Zds.'+vord+'.mseed',format='MSEED')
    return Nss,Ess,Zss,Nds,Eds,Zds


def list2stream(st_list):
    
    from obspy import Stream
    
    st=Stream()
    for k in range(len(st_list)):
        st += st_list[k]
    return st
    


def load_fakequakes_tele_synthetics(home,project_name,fault_name,model_name,GF_list_teleseismic,G_from_file,G_name,decimation_factor):
    ''''
    Load the miniseed files with all the synthetics
    '''
    from numpy import genfromtxt,loadtxt
    from obspy import read,Stream

    if G_from_file==True: #load from file
        Eds=read(home+project_name+'/GFs/matrices/'+G_name+'.Eds.tele.mseed')
        Nds=read(home+project_name+'/GFs/matrices/'+G_name+'.Nds.tele.mseed')
        Zds=read(home+project_name+'/GFs/matrices/'+G_name+'.Zds.tele.mseed')
        Ess=read(home+project_name+'/GFs/matrices/'+G_name+'.Ess.tele.mseed')
        Nss=read(home+project_name+'/GFs/matrices/'+G_name+'.Nss.tele.mseed')
        Zss=read(home+project_name+'/GFs/matrices/'+G_name+'.Zss.tele.mseed')
    else: #assemble G one data type at a time, just displacememnt right now
        #Load station info
        station_file=home+project_name+'/data/station_info/'+GF_list_teleseismic
        staname=genfromtxt(station_file,dtype="U",usecols=0)
        Nsta=len(staname)
        #Load fault model
        source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
        Nfaults=source.shape[0] #Number of subfaults
        kindex=0
        
        for ksta in range(Nsta):
            
            print('Reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
            
            for kfault in range(Nfaults):
                
                #Get subfault GF directory
                nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                strdepth='%.4f' % source[kfault,3]
                syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                
                #Get synthetics
                if kfault==0 and ksta==0: #It's the first one, initalize stream object
                    SS=read(syn_path+staname[ksta]+'.'+nfault+'.SS.mseed')
                    DS=read(syn_path+staname[ksta]+'.'+nfault+'.DS.mseed')
                    
                    for trace in SS:
                        
                        if decimation_factor >1:
                            trace.decimate(factor=decimation_factor)    
                        
                        
                        if trace.stats.channel=='BXE':
                            Ess=Stream(trace)
                        if trace.stats.channel=='BXN':
                            Nss=Stream(trace)
                        if trace.stats.channel=='BXZ':
                            Zss=Stream(trace)
                        
                    for trace in DS:
                        
                        if decimation_factor >1:
                            trace.decimate(factor=decimation_factor)    
                        
                        if trace.stats.channel=='BXE':
                            Eds=Stream(trace)
                        if trace.stats.channel=='BXN':
                            Nds=Stream(trace)
                        if trace.stats.channel=='BXZ':
                            Zds=Stream(trace)
                            

                else: #Just add to stream object
                   
                    SS=read(syn_path+staname[ksta]+'.'+nfault+'.SS.mseed')
                    DS=read(syn_path+staname[ksta]+'.'+nfault+'.DS.mseed')
                    
                    for trace in SS:
                        
                        if decimation_factor >1:
                            trace.decimate(factor=decimation_factor)  
                        
                        if trace.stats.channel=='BXE':
                            Ess+=trace
                        if trace.stats.channel=='BXN':
                            Nss+=trace
                        if trace.stats.channel=='BXZ':
                            Zss+=trace
                        
                    for trace in DS:
                        
                        if decimation_factor >1:
                            trace.decimate(factor=decimation_factor)  
                        
                        if trace.stats.channel=='BXE':
                            Eds+=trace
                        if trace.stats.channel=='BXN':
                            Nds+=trace
                        if trace.stats.channel=='BXZ':
                            Zds+=trace
                kindex+=1
        print('Writting synthetics to miniSEED, hang on this might take a minute or two.')
        Ess.write(home+project_name+'/GFs/matrices/'+G_name+'.Ess.tele.mseed',format='MSEED')
        Nss.write(home+project_name+'/GFs/matrices/'+G_name+'.Nss.tele.mseed',format='MSEED')
        Zss.write(home+project_name+'/GFs/matrices/'+G_name+'.Zss.tele.mseed',format='MSEED')
        Eds.write(home+project_name+'/GFs/matrices/'+G_name+'.Eds.tele.mseed',format='MSEED')
        Nds.write(home+project_name+'/GFs/matrices/'+G_name+'.Nds.tele.mseed',format='MSEED')
        Zds.write(home+project_name+'/GFs/matrices/'+G_name+'.Zds.tele.mseed',format='MSEED')
    return Nss,Ess,Zss,Nds,Eds,Zds





def get_fakequakes_G_and_m(Gimpulse,home,project_name,rupture_name,time_epi,GF_list,epicenter,NFFT,
                source_time_function,stf_falloff_rate,dt=1.0,zeta=0.2,forward=False,reconvolution=False,
                old_stf='prem_i_2s',ncpus=1):
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
    
    from numpy import genfromtxt,convolve,where,zeros,arange,unique,r_,sort,ones
    from numpy import expand_dims,squeeze,roll,float64,array
    import dask.array as da

    if forward==True:
        source=genfromtxt(home+project_name+'/forward_models/'+rupture_name)
    else:
        source=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)
    rise_times=source[:,7]
    fraction=source[:,6]
    rupture_onset=source[:,12]
    
    # Gaussian STF hard coded values    WARNING!!!!
    old_rise_time = 2.1 #This is the prem_i_2s hard coded value
    time_offset_gauss=100 #Hard coded for now. Long enough for any conceivable rise time
    #####################
    
    
    
    #How many unique faults?
    Nfaults=len(unique(source[:,0]))
    
    #How many subfaults are non-zero?
    i_non_zero_original=where(rise_times>0)[0]
    N_non_zero=len(i_non_zero_original)
    rise_times = rise_times[i_non_zero_original]
    fraction = fraction[i_non_zero_original]
    rupture_onset = rupture_onset[i_non_zero_original]
    
    #Convert rupture onsets to itneger number of samples
    rupture_onset = (rupture_onset/dt).astype('int')
    
    #Now convert to indices according to how G matrix is ordered
    i_non_zero = sort(r_[i_non_zero_original*2,i_non_zero_original*2+1])
    
    #Stations
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    Nsta=len(staname)
    
    #Extract only impulses from non-zero fautls and put into
    #Dask array Convert impulses matrix to dask array
    Gimpulse = da.from_array(Gimpulse[:,i_non_zero],chunks = (NFFT*Nsta*3,2))
    
    #Make array of slip rate functions, it will be one matrix for a gingle sta 
    #and all subfaults simualtenously
    slip_rates = ones((NFFT,N_non_zero))
    
    total_time = (NFFT-1) * dt
    for ksub in range(N_non_zero):
        
        tau_r = rise_times[ksub]
        ji_fraction = fraction[ksub] 
        
        if reconvolution == False:
            sr = build_source_time_function(tau_r,dt,total_time,stf_type=source_time_function,
                            zeta=zeta,dreger_falloff_rate=stf_falloff_rate,scale=True,scale_value=dt,
                            ji_fraction = ji_fraction)
        else:
            sr = build_source_time_function(tau_r,dt,total_time,stf_type='gauss_prem_i_2s',
                            time_offset_gauss=time_offset_gauss,scale=True,scale_value=1.0,quiet=True)
               
        t_stf,slip_rates[:,ksub] = expand_dims(sr,1)
        
    #COnvert slip rates and onsets to dask array
    rupture_onset = da.from_array(rupture_onset,chunks = (1,))
    slip_rates = da.from_array(slip_rates,chunks = (NFFT,1))
    
    
    #Are we reconvovling?
    if reconvolution == True:
        #Need old STF
        t_stf,old_stf=build_source_time_function(old_rise_time,dt,total_time,
                stf_type='gauss_prem_i_2s',time_offset_gauss=time_offset_gauss,
                scale=True,scale_value=1.0)
    else:
        old_stf = None

    print('... ... done calculating slip rate functions, block convolution is next')
    
    
    #Define block convolution and trim function to be mapped with DASK
    def block_convolution_and_shift(gf_block,slip_rate_block,onset_block,NFFT,
                    old_slip_rate_function = None,reconvolution = False):
        '''
        gf_block is the block from the impulse response GFs
        slip_rate_block  is the block from the slip rate functions
        '''
        
        start_row = 0 #These two track which bit of the array to extract
        end_row = NFFT
        
        
        #Define slip rate functions
        slip_rate_function = squeeze(slip_rate_block) #New STF
        if reconvolution == True:
             old_slip_rate_function = squeeze(old_slip_rate_function)
        
        Gout = ones(gf_block.shape)
        
        for ksta in range(Nsta):
            for kcomponent in range(3): # Loops over N,E,Z
                for krake in range(2): #loops over ss and ds
            
                    gf = gf_block[start_row:end_row,krake]
                    
                    #remove extraneous dimensions
                    gf = squeeze(gf)
                   
                    
                    #Roll the gf forward by rupture_onset samples
                    gf = roll(gf,onset_block)
                    #Zero out the first rolled samples
                    gf[0:onset_block[0]] = 0 #The sonet block is really jsut a single number
                    

                    
                    #Convovle impulse resp and slip rate function
                    if reconvolution == False :  #COnvovle new STF with impulse
                        gf = convolve(gf,slip_rate_function)[0:NFFT]
                    else:  #Reconvolution, need both new STF and old STF
                        gf = gauss_reconvolution(gf,slip_rate_function,old_slip_rate_function)
                        
        
                    #save for output
                    Gout[start_row:end_row,krake] = gf
                
                #Update counters
                start_row += NFFT
                end_row += NFFT   
        
        return Gout
    
    G = da.map_blocks(block_convolution_and_shift, Gimpulse, slip_rates,
            rupture_onset,NFFT,old_stf,reconvolution,chunks=(NFFT*Nsta*3,2),
            dtype=float64).persist(num_workers=ncpus)

    
    #Get slip model vector
    m=zeros((N_non_zero*2,1))
    iss=arange(0,len(m),2)
    ids=arange(1,len(m),2)
    m[iss,0]=source[i_non_zero_original,8]
    m[ids,0]=source[i_non_zero_original,9]

    return m,array(G)




def get_fakequakes_G_and_m_old_pre_DASK(Gimpulse,home,project_name,rupture_name,time_epi,GF_list,epicenter,NFFT,
                source_time_function,stf_falloff_rate,zeta=0.2,forward=False,reconvolution=False,old_stf='prem_i_2s'):
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
    
    from numpy import genfromtxt,convolve,where,zeros,arange,unique
    import gc
    import dask.array as da


    if forward==True:
        source=genfromtxt(home+project_name+'/forward_models/'+rupture_name)
    else:
        source=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)
    rise_times=source[:,7]
    rupture_onset=source[:,12]
    
    #How many unique faults?
    Nfaults=len(unique(source[:,0]))
    
    #How many subfaults are non-zero?
    i_non_zero=where(rise_times>0)[0]
    N_non_zero=len(i_non_zero)
    
    #Stations
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    Nsta=len(staname)
    
    #Initalize G matrix
    G=zeros((NFFT*3*Nsta,N_non_zero*2))
    
    #Place synthetics in correct place in Gmatrix
    matrix_pos=0 #tracks where in matrix synths are placed
    read_start=0  #Which trace to start reading from
    for ksta in range(Nsta):
        print('... ... working on station '+str(ksta+1)+' of '+str(Nsta))
        
        for ksource in range(len(i_non_zero)):

            if ksource % 50 == 0:
                print('... ... ... subfault %d of %d' % (ksource,len(i_non_zero)))
            
            #Get synthetics
            nss=Nss[read_start+i_non_zero[ksource]].copy()
            ess=Ess[read_start+i_non_zero[ksource]].copy()
            zss=Zss[read_start+i_non_zero[ksource]].copy()
            nds=Nds[read_start+i_non_zero[ksource]].copy()
            eds=Eds[read_start+i_non_zero[ksource]].copy()
            zds=Zds[read_start+i_non_zero[ksource]].copy()
            #Delay synthetics by rupture onset
            tdelay=rupture_onset[i_non_zero[ksource]]
            nss,ess,zss,nds,eds,zds=tshift_trace(nss,ess,zss,nds,eds,zds,tdelay,time_epi,NFFT)
            #Convolve with source time function
            rise=rise_times[i_non_zero[ksource]]
            #Make sure rise time is a multiple of dt
            dt=nss.stats.delta
            rise=round(rise/dt)*nss.stats.delta
            if rise<(2.*dt): #Otherwise get nan's in STF
                rise=2.*dt
            
            total_time=NFFT*dt-dt
            
            if reconvolution==False: #Do a straight up convolution with new slip rate function
                
                #Build the new slip rate function
                t_stf,stf=build_source_time_function(rise,dt,total_time,stf_type=source_time_function,zeta=zeta,dreger_falloff_rate=stf_falloff_rate)

                nss.data=convolve(nss.data,stf)[0:NFFT]
                ess.data=convolve(ess.data,stf)[0:NFFT]
                zss.data=convolve(zss.data,stf)[0:NFFT]
                nds.data=convolve(nds.data,stf)[0:NFFT]
                eds.data=convolve(eds.data,stf)[0:NFFT]
                zds.data=convolve(zds.data,stf)[0:NFFT]
            else: #reconvovleby using old stf
            

                
                #Need old STF
                old_rise_time = 2.1 #This is the prem_i_2s hard coded value
                time_offset_gauss=100 #Hard coded for now. Long enough for any conceivable rise time
                t_stf,old_stf=build_source_time_function(old_rise_time,dt,total_time,stf_type='gauss_prem_i_2s',time_offset_gauss=time_offset_gauss,scale=True,scale_value=1.0)
                
                #Now the new STF
                t_stf,new_stf=build_source_time_function(rise,dt,total_time,stf_type='gauss_prem_i_2s',time_offset_gauss=time_offset_gauss,scale=True,scale_value=1.0,quiet=True)
                
                #Reconvolutions
                nss.data=gauss_reconvolution(nss.data,new_stf,old_stf)
                ess.data=gauss_reconvolution(ess.data,new_stf,old_stf)
                zss.data=gauss_reconvolution(zss.data,new_stf,old_stf)
                nds.data=gauss_reconvolution(nds.data,new_stf,old_stf)
                eds.data=gauss_reconvolution(eds.data,new_stf,old_stf)
                zds.data=gauss_reconvolution(zds.data,new_stf,old_stf)
                
                
            #Place in matrix
            G[matrix_pos:matrix_pos+NFFT,2*ksource]=nss.data
            G[matrix_pos:matrix_pos+NFFT,2*ksource+1]=nds.data
            G[matrix_pos+NFFT:matrix_pos+2*NFFT,2*ksource]=ess.data
            G[matrix_pos+NFFT:matrix_pos+2*NFFT,2*ksource+1]=eds.data
            G[matrix_pos+2*NFFT:matrix_pos+3*NFFT,2*ksource]=zss.data
            G[matrix_pos+2*NFFT:matrix_pos+3*NFFT,2*ksource+1]=zds.data
            
            #Having some sort of memory issue delete streams and garbage collect to avoid it
            del nss,ess,zss,nds,eds,zds
            gc.collect()
            
        matrix_pos+=3*NFFT
        read_start+=Nfaults
        #Get slip model vector
    m=zeros((N_non_zero*2,1))
    iss=arange(0,len(m),2)
    ids=arange(1,len(m),2)
    m[iss,0]=source[i_non_zero,8]
    m[ids,0]=source[i_non_zero,9]

    return m,G








def gauss_reconvolution(data,new_stf,old_stf):
    
    '''
    Used if original slip rate function was NOT  an impulse, then you need to "reconvolve"
    as per the SRL paper on Syngine/Instaseis (kirscher et al., 2017?). THis currently
    will only work for teleseismics and the slip rate functions used by instaseis
    which are all Gaussians
    '''
    
    from numpy import fft,real
    
    #fft STFs and data
    Data=fft.fft(data) 
    STF_new=fft.fft(new_stf)
    STF_old=fft.fft(old_stf)
    
    #New seismogram
    Data_new=(Data*STF_new)/STF_old
    data_new=real(fft.ifft(Data_new))
    
    return data_new


def get_fakequakes_G_and_m_dynGF(Nss,Ess,Zss,Nds,Eds,Zds,home,project_name,rupture_name,time_epi,STA,GF_list,epicenter,NFFT,
                source_time_function,stf_falloff_rate,forward=False):
    '''
    Assemble Green functions matrix. If requested will parse all available synthetics on file and build the matrix.
    Otherwise, if it exists, it will be loaded from file 
    
    IN:
        home: Home directory location
        project_name: Name of the problem
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        STA: Dictionary of station's index
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

    from numpy import genfromtxt,loadtxt,convolve,where,zeros,arange,unique,save
    import numpy as np

    if forward==True:
        source=genfromtxt(home+project_name+'/forward_models/'+rupture_name)
    else:
        source=genfromtxt(home+project_name+'/output/ruptures/'+rupture_name)
    rise_times=source[:,7]
    rupture_onset=source[:,12]

    #How many unique faults?
    Nfaults=len(unique(source[:,0]))

    #How many subfaults are non-zero?
    i_non_zero=where(rise_times>0)[0]
    N_non_zero=len(i_non_zero)

    #Stations
    station_file=home+project_name+'/data/station_info/'+GF_list
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    staname=[n_name.decode() for n_name in staname]
    Nsta=len(staname)

    #Initalize G matrix
    G=zeros((NFFT*3*Nsta,N_non_zero*2))

    #Place synthetics in correct place in Gmatrix
    matrix_pos=0 #tracks where in matrix synths are placed
    #read_start=0  #Which trace to start reading from
    for ksta in range(Nsta):
        print('... ... working on station '+str(ksta)+' of '+str(Nsta))
        read_start=STA[staname[ksta]]*Nfaults #Which trace to start reading from. Depending on the original GF_list order
        for ksource in range(len(i_non_zero)):

            #Get synthetics
            nss=Nss[read_start+i_non_zero[ksource]].copy()
            ess=Ess[read_start+i_non_zero[ksource]].copy()
            zss=Zss[read_start+i_non_zero[ksource]].copy()
            nds=Nds[read_start+i_non_zero[ksource]].copy()
            eds=Eds[read_start+i_non_zero[ksource]].copy()
            zds=Zds[read_start+i_non_zero[ksource]].copy()
            #Delay synthetics by rupture onset
            tdelay=rupture_onset[i_non_zero[ksource]]
            nss,ess,zss,nds,eds,zds=tshift_trace(nss,ess,zss,nds,eds,zds,tdelay,time_epi,NFFT)
            #Convolve with source time function
            rise=rise_times[i_non_zero[ksource]]
            #Make sure rise time is a multiple of dt
            dt=nss.stats.delta
            rise=round(rise/dt)*nss.stats.delta
            if rise<(2.*dt): #Otherwise get nan's in STF
                rise=2.*dt
            total_time=NFFT*dt
            t_stf,stf=build_source_time_function(rise,dt,total_time,stf_type=source_time_function,dreger_falloff_rate=stf_falloff_rate)

            nss.data=convolve(nss.data,stf)[0:NFFT]
            ess.data=convolve(ess.data,stf)[0:NFFT]
            zss.data=convolve(zss.data,stf)[0:NFFT]
            nds.data=convolve(nds.data,stf)[0:NFFT]
            eds.data=convolve(eds.data,stf)[0:NFFT]
            zds.data=convolve(zds.data,stf)[0:NFFT]
            #Place in matrix
            G[matrix_pos:matrix_pos+NFFT,2*ksource]=nss.data
            G[matrix_pos:matrix_pos+NFFT,2*ksource+1]=nds.data
            G[matrix_pos+NFFT:matrix_pos+2*NFFT,2*ksource]=ess.data
            G[matrix_pos+NFFT:matrix_pos+2*NFFT,2*ksource+1]=eds.data
            G[matrix_pos+2*NFFT:matrix_pos+3*NFFT,2*ksource]=zss.data
            G[matrix_pos+2*NFFT:matrix_pos+3*NFFT,2*ksource+1]=zds.data

        matrix_pos+=3*NFFT
        #read_start+=Nfaults  #This is not always true, unless the GF_list is in the same order and use the same stations
        #Get slip model vector
    m=zeros((N_non_zero*2,1))
    iss=arange(0,len(m),2)
    ids=arange(1,len(m),2)
    m[iss,0]=source[i_non_zero,8]
    m[ids,0]=source[i_non_zero,9]

    return m,G





def write_fakequakes_waveforms(home,project_name,rupture_name,waveforms,GF_list,NFFT,time_epi,dt):
    '''
    write output from fakequakes run
    '''
    
    from os import path,makedirs
    from numpy import genfromtxt,squeeze,sqrt
    from obspy import Stream,Trace
    from pyproj import Geod
    
    #Where am I writting to?
    #old
#    rupture=rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new
    rupture = rupture_name.rsplit('.', 1)[0]
    
    directory=home+project_name+'/output/waveforms/'+rupture+'/'
    
    #Check if dir exists if not then create it
    if not path.exists(directory):
        makedirs(directory)
    
    # Get station list and coords
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U') 
    lon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=1)  
    lat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=2) 
    # Get rupture information to store in waveform SAC files
    epicenter,time_epi,Mw=read_fakequakes_hypo_time(home,project_name,rupture_name,get_Mw=True)
    
    #Projection object for distance calculations
    g=Geod(ellps='WGS84')
    
    line_out='# Station, lon, lat, N[m], E[m], Up[m], PGD[m]\n'
    #Parse waveforms vector
    read_pos=0
    for ksta in range(len(sta)):
        # Init waveforms and metadata
        n=Stream(Trace())
        n[0].stats.starttime=time_epi
        n[0].stats.delta=dt
        n[0].stats.station=sta[ksta]
        az,backaz,dist_in_m=g.inv(epicenter[0],epicenter[1],lon[ksta],lat[ksta])
        dist_in_km=(1./1000)*dist_in_m
        n[0].stats.update({'sac':{'stlo':lon[ksta],'stla':lat[ksta],'evlo':epicenter[0],'evla':epicenter[1],'evdp':epicenter[2],'dist':dist_in_km,'az':az,'baz':backaz,'mag':Mw}})
#        n[0].stats.latitude=lat[ksta]  #these don't appear to save properly in SAC format anyway
#        n[0].stats.longitude=lon[ksta]
        # Copy to the other components
        e=n.copy()
        z=n.copy()
        # Get the actual data
        n[0].data=squeeze(waveforms[read_pos:read_pos+NFFT])
        e[0].data=squeeze(waveforms[read_pos+NFFT:read_pos+2*NFFT])
        z[0].data=squeeze(waveforms[read_pos+2*NFFT:read_pos+3*NFFT])
        # Extract coseismic offsets
        n_offset=n[0].data[-10:].mean()
        e_offset=e[0].data[-10:].mean()
        z_offset=z[0].data[-10:].mean()
        # Get PGD
        pgd=sqrt(n[0].data**2+e[0].data**2+z[0].data**2).max()
        # Summary file
        line_out+='%s\t%10.4f\t%10.4f\t%.5e\t%.6e\t%.5e\t%.5e\n' % (sta[ksta],lon[ksta],lat[ksta],n_offset,e_offset,z_offset,pgd)
        # write to file
        n.write(directory+sta[ksta]+'.LYN.sac',format='SAC')
        e.write(directory+sta[ksta]+'.LYE.sac',format='SAC')
        z.write(directory+sta[ksta]+'.LYZ.sac',format='SAC')
        # update counter
        read_pos+=3*NFFT
    f=open(directory+'_summary.'+rupture+'.txt','w')
    f.write(line_out)
    f.close()
        
        
       


def coseismics(home,project_name,rupture_name,station_file,hot_start=None):
    '''
    This routine will take synthetics and apply a static slip dsitibution. It will 
    linearly superimpose the synthetic coseismic from each subfault. Output will be
    a single ascii file witht he 3 coseismic offsets (NEU) for each station defined 
    in the station_file. Depending on the specified rake angle at each subfault the 
    code will compute the contribution to dip and strike slip directions. It will 
    also compute the moment at that subfault and scale it according to the unit 
    amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
       
    OUT:
        Nothing
    '''
    from numpy import loadtxt,genfromtxt,array,savetxt,unique,where
    
    print('Solving for static problem')
    #Output where?
    outpath=home+project_name+'/output/forward_models/'
    #load source
    source=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    #Load stations
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    #Get unique sources
    source_id=unique(source[:,0])
    #Loop over stations
    if hot_start is None:
        hot_start=0
    for ksta in range(hot_start,len(staname)):
        #Initalize output
        n=array([0])
        e=array([0])
        z=array([0])
        sta=staname[ksta]
        print('Working on station '+staname[ksta]+' ('+str(ksta+1)+'/'+str(len(staname))+')')
        #Loop over sources
        for k in range(len(source_id)):
            print(k)
            #Get subfault parameters
            nfault='subfault'+str(int(source_id[k])).rjust(4,'0')
            ifault=where(source[:,0]==source_id[k])[0]
            ss_slip=source[ifault,8].sum()
            ds_slip=source[ifault,9].sum()
            print('ds_slip='+str(ds_slip))
            #Where's the data
            syn_path=home+project_name+'/GFs/static/'
            #Get synthetics
            try:
                coseis_ss=loadtxt(syn_path+sta+'.'+nfault+'.SS.static.neu')
            except:
                coseis_ss=array([0,0,0])
            nss=coseis_ss[0]
            ess=coseis_ss[1]
            zss=coseis_ss[2]
            try:
                coseis_ds=loadtxt(syn_path+sta+'.'+nfault+'.DS.static.neu')
            except:
                coseis_ds=array([0,0,0])
            nds=coseis_ds[0]
            eds=coseis_ds[1]
            zds=coseis_ds[2]
            print('zds='+str(zds))
            #get rake contribution and moment multiplier
            etotal=ds_slip*eds+ss_slip*ess
            ntotal=ds_slip*nds+ss_slip*nss
            ztotal=ds_slip*zds+ss_slip*zss
            print('ztotal='+str(ztotal))
            #Add to previous subfault's results
            e=e+etotal
            n=n+ntotal
            z=z+ztotal
            print('n='+str(n))
            print('e='+str(e))
            print('z='+str(z))
        #Save results
        savetxt(outpath+sta+'.static.neu',(n,e,z))


def coseismics_matrix(home,project_name,rupture_name,station_file,G_from_file,
                      model_name,G_name,return_G=False,G=None):
    '''
    This routine will take synthetics and apply a static slip dsitibution. It will 
    linearly superimpose the synthetic coseismic from each subfault. Output will be
    a single ascii file witht he 3 coseismic offsets (NEU) for each station defined 
    in the station_file. Depending on the specified rake angle at each subfault the 
    code will compute the contribution to dip and strike slip directions. It will 
    also compute the moment at that subfault and scale it according to the unit 
    amount of momeent (1e15 N-m)
    
    IN:
        home: Home directory
        project_name: Name of the problem
        rupture_name: Name of rupture description file
        station_file: File with coordinates of stations
        model_Name: Name of Earth structure model file
       
    OUT:
        Nothing
    '''
    from numpy import loadtxt,genfromtxt,array,savetxt,unique,where,zeros,load,save,arange
    import os
    from glob import glob
    
    print('... solving for static problem')
    #Output where?
    outpath=home+project_name+'/output/statics/'+rupture_name.replace('.rupt','')+'/'
    
    #check if output folder exists, if not make it
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    #load source
    source=loadtxt(home+project_name+'/output/ruptures/'+rupture_name,ndmin=2)
    print('... ... source geometry is '+home+project_name+'/output/ruptures/'+rupture_name)
    
    #Load stations
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    sta_lonlat=genfromtxt(station_file,usecols=[1,2])
    
    #Get unique sources
    source_id=unique(source[:,0])
    m=zeros((2*len(source_id),1))
    #Assemble source parameters vector
    for ksource in range(len(source_id)):

        ifault=where(source[:,0]==source_id[ksource])[0]
        #Combine into model vector
        ss_slip=source[ifault,8].sum()
        ds_slip=source[ifault,9].sum()
        m[2*ksource]=ss_slip
        m[2*ksource+1]=ds_slip
    
    
    #Do I need to make G or can I load it from file??
    G_name=home+project_name+'/GFs/matrices/'+G_name
    
    if G_from_file==True and G is None: #load from file
        
        if G_name[-3:]!='npy':
            G_name=G_name+'.npy'
        print('... .... loading G from file '+G_name)
        G=load(G_name)
    
    elif G_from_file==True and G is not None:
        print('... ... G already provided, do not reload ...')
    
    elif G_from_file==False and G is None:
        
        #initalize matrices
        G=zeros((3*len(staname),2*len(source_id)))
        print('... ... Assembling GFs matrix from scratch')
        
            
        for ksource in range(len(source_id)):
            
            if ksource % 50 == 0:
                print('... ... Loading source %i of %i' %(ksource,len(source_id)))
            
            #Get subfault parameters
            nfault='sub'+str(int(source_id[ksource])).rjust(4,'0')
            nfault2='subfault'+str(int(source_id[ksource])).rjust(4,'0')
            #Where's the synthetic data
            syn_path=glob(home+project_name+'/GFs/static/'+model_name+'*'+nfault)[0]
            
            #Get synthetics
            try:
                coseis_ss=genfromtxt(syn_path+'/'+nfault2+'.SS.static.neu',usecols=[1,2,3])
            except:
                coseis_ss=array([0,0,0])
            
            Nsites = len(coseis_ss)
            
            nss=coseis_ss[:,0]
            ess=coseis_ss[:,1]
            zss=coseis_ss[:,2]
            
            try:
                coseis_ds=genfromtxt(syn_path+'/'+nfault2+'.DS.static.neu',usecols=[1,2,3])
            except:
                coseis_ds=array([0,0,0])
            
            nds=coseis_ds[:,0]
            eds=coseis_ds[:,1]
            zds=coseis_ds[:,2]
            
            #Put in matrix
            inorth=arange(0,Nsites*3,3)
            ieast=arange(1,Nsites*3,3)
            iup=arange(2,Nsites*3,3)
            
            #North
            G[inorth,2*ksource]=nss  ;  G[inorth,2*ksource+1]=nds
            #East
            G[ieast,2*ksource]=ess  ;  G[ieast,2*ksource+1]=eds
            #East
            G[iup,2*ksource]=zss  ;  G[iup,2*ksource+1]=zds
        
        #Save G matrix
        print('Saving GF matrix to '+G_name+' this might take just a second...')
        save(G_name,G)
    
    #Now go on to matrix multiply and save solutions
    print('... matrix multiplying and saving output...DONE')
    d=G.dot(m)
    
    #write to file
    out_file = outpath+'statics.neu'
    f=open(out_file,'w')
    f.write('# sta,lon,lat,n(m),e(m),z(m)\n')
            
    for ksta in range(len(staname)):
        sta=staname[ksta]
        n=d[3*ksta]
        e=d[3*ksta+1]
        z=d[3*ksta+2]
        
        f.write('%s\t%.4f\t%.4f\t%.4e\t%.4e\t%.4e\n' % (sta,sta_lonlat[ksta,0],sta_lonlat[ksta,1],n,e,z))
    f.close()
    
    #return G?
    if return_G==True:
        return G
    

def coseismics_fakequakes(home,project_name,GF_list,master_G_from_file,G_name,
                          model_name,rupture_list):
    '''
    Make static offsets for all fakequakes ruptures
    '''

    from numpy import genfromtxt,where
    
    #load rupture list
    all_sources=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='U')
    
    #make temp GF_list with onlys tations that are flagged as statics
    stations=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U')
    station_lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    statics_flag=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=3)
    
    #Keep only relevant sites
    i=where(statics_flag==1)[0]
    stations=stations[i]
    station_lonlat=station_lonlat[i]
    
    station_file=home+project_name+'/data/station_info/coseismics.tmp'
    f=open(station_file,'w')
    for ksta in range(len(stations)):
        f.write('%s\t%.4f\t%.4f\n' % (stations[ksta],station_lonlat[ksta,0],station_lonlat[ksta,1]))
    f.close()
        
    #loop over sources
    for krupt in range(len(all_sources)):
        
        #Run coseismcis on one rupture
        rupture_name=all_sources[krupt]
        print('Working on rupture '+rupture_name)
        
        if master_G_from_file==True:
            if krupt==0:
                G_from_file=True
                G=coseismics_matrix(home,project_name,rupture_name,'coseismics.tmp',G_from_file,
                                      model_name,G_name,return_G=True,G=None)
            else:
                G_from_file=True
                coseismics_matrix(home,project_name,rupture_name,'coseismics.tmp',G_from_file,
                                  model_name,G_name,return_G=False,G=G)
        else:
        #Only make G matrix the first time, otehrwise load it from file
            if krupt==0:
                G_from_file=False
                G=coseismics_matrix(home,project_name,rupture_name,'coseismics.tmp',G_from_file,
                                    model_name,G_name,return_G=True,G=None)
            else:
                G_from_file=True
                coseismics_matrix(home,project_name,rupture_name,'coseismics.tmp',G_from_file,
                                  model_name,G_name,return_G=False,G=G)
            

        





def tsunami_waveforms(home,project_name,fault_name,rupture_name,station_file,model_name,run_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,coord_type,decimate,lowpass,resample,beta):            
    '''
    Forward compute tsunami waveforms the right way, load the GF matrix and just multiply by fault model
    '''
    from mudpy.inverse import getG
    from numpy import genfromtxt,zeros,arange,deg2rad,cos,sin,vstack,array,where
    from obspy import read

    #Read GFs
    G=getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,coord_type,decimate,lowpass)
    #Read rupture model and convert into vector
    ss=genfromtxt(home+project_name+'/forward_models/'+rupture_name,usecols=8)
    ds=genfromtxt(home+project_name+'/forward_models/'+rupture_name,usecols=9)
    #Rotate by beta
    if beta != None:
        beta_rot=deg2rad(beta)
        R=array([[cos(beta_rot),sin(beta_rot)],[-sin(beta_rot),cos(beta_rot)]])
        rot=R.dot(vstack((ss,ds)))
        ss_rot=rot[0,:]
        ds_rot=rot[1,:]
    #Assemble into column vector
    m=zeros(len(ss)*2)
    iss=arange(0,len(m),2)
    ids=arange(1,len(m),2)
    m[iss]=ss_rot
    m[ids]=ds_rot
    #Multiply
    dtsun=G.dot(m)
    #Write to file (Modified from inverse.write_synthetics)
    #Read gf file and decide what needs to get loaded
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=[0],skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[3,4,5,6,7],skip_header=1,dtype='f8')
    GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],skip_header=1,dtype='U')
    #Separate into its constituent parts (statics,displacaments, velocities, etc...)
    kinsert=0
    kgf=3
    i=where(GF[:,kgf]==1)[0]
    if len(i)>0:
        for ksta in range(len(i)):
            sta=stations[i[ksta]]
            tsun=read(GFfiles[i[ksta],kgf])
            npts=tsun[0].stats.npts
            synth=tsun.copy()
            synth[0].data=dtsun[kinsert:kinsert+npts]
            kinsert+=npts
            synth.write(home+project_name+'/output/forward_models/'+run_name+'.'+sta+'.tsun',format='SAC')

            
                        
                                 
                                                            
                        
                                                
def move_seafloor(home,project_name,run_name,topo_dx_file,topo_dy_file,tgf_file,
        outname,time_epi,tsun_dt,maxt,coast_file,topo_effect=True,variance=None,static=False):
    '''
    Create moving topography input files for geoclaw
    '''
    import datetime
    from numpy import genfromtxt,zeros,meshgrid,ones,c_,savetxt,delete,where,nan,argmin,arange
    from obspy import read
    from netCDF4 import Dataset
    from scipy.interpolate import griddata
    from mudpy.inverse import interp_and_resample,grd2xyz
    from scipy.ndimage.filters import gaussian_filter


    #Get station names
    sta=genfromtxt(home+project_name+'/data/station_info/'+tgf_file)
    stanames=genfromtxt(home+project_name+'/data/station_info/'+tgf_file,usecols=0,dtype='U')
    #lon=360+sta[:,1]
    lon=sta[:,1]
    print('correcting longitude')
    print(lon[0])
    lat=sta[:,2]
    #Get fault file
    #f=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #Where is the data
    data_dir=home+project_name+'/output/forward_models/'
    #Define time deltas
    td_max=datetime.timedelta(seconds=maxt)
    #Maximum tiem to be modeled
    tmax=time_epi+td_max
    #Read derivatives
    bathy_dx=Dataset(topo_dx_file,'r',format='NETCDF4')
    zdx=bathy_dx.variables['z'][:].data
    bathy_dy=Dataset(topo_dy_file,'r',format='NETCDF4')
    zdy=bathy_dy.variables['z'][:].data
    #Make interpolation quantities
    try:
        loni=bathy_dx.variables['x'][:]
        lati=bathy_dx.variables['y'][:]
    except:
        loni=bathy_dx.variables['lon'][:]
        lati=bathy_dx.variables['lat'][:]
    loni,lati=meshgrid(loni,lati)
    #Read slope file
    kwrite=0
    idelete=[]
    for ksta in range(len(sta)):
        if ksta%500==0:
            print('... ... working on seafloor grid point '+str(ksta)+' of '+str(len(sta)))
        try: #If no data then delete
            if static==False: #We're reading waveforms
                #e=read(data_dir+run_name+'.'+rjust(str(int(sta[ksta,0])),4,'0')+'.disp.e')
                #n=read(data_dir+run_name+'.'+rjust(str(int(sta[ksta,0])),4,'0')+'.disp.n')
                #u=read(data_dir+run_name+'.'+rjust(str(int(sta[ksta,0])),4,'0')+'.disp.z')
                #e=read(data_dir+run_name+'.'+stanames[ksta]+'.disp.e')
                #n=read(data_dir+run_name+'.'+stanames[ksta]+'.disp.n')
                #u=read(data_dir+run_name+'.'+stanames[ksta]+'.disp.u')
                e=read(data_dir+stanames[ksta]+'.disp.e')
                n=read(data_dir+stanames[ksta]+'.disp.n')
                u=read(data_dir+stanames[ksta]+'.disp.u')
                e=interp_and_resample(e,1.0,time_epi)
                n=interp_and_resample(n,1.0,time_epi)
                u=interp_and_resample(u,1.0,time_epi)
                #Keep only data between time_epi and tmax
                e.trim(time_epi,tmax,fill_value=e[0].data[-1],pad=True)
                n.trim(time_epi,tmax,fill_value=n[0].data[-1],pad=True)
                u.trim(time_epi,tmax,fill_value=u[0].data[-1],pad=True)
                #Initalize matrices
                if ksta==0:
                    emat=zeros((n[0].stats.npts,len(sta)))
                    nmat=emat.copy()
                    umat=emat.copy()
                #Populate matrix
                emat[:,kwrite]=e[0].data
                nmat[:,kwrite]=n[0].data
                umat[:,kwrite]=u[0].data
            else:
                neu=genfromtxt(data_dir+stanames[ksta]+'.static.neu')
                n=neu[0]
                e=neu[1]
                u=neu[2]
                tsun_dt=1.0
                maxt=1.0
                if ksta==0:
                    emat=zeros((1,len(sta)))
                    nmat=emat.copy()
                    umat=emat.copy()
                #Populate matrix
                emat[:,kwrite]=e
                nmat[:,kwrite]=n
                umat[:,kwrite]=u
            kwrite+=1
        except: #Data was missing, delete from lat,lon
            pass
            print('No data for station '+str(ksta)+', deleting from coordinates list')
            idelete.append(ksta)
    #Clean up missing data
    if len(idelete)!=0:
        lat=delete(lat,idelete)
        lon=delete(lon,idelete)
        emat=emat[:,:-len(idelete)]
        nmat=nmat[:,:-len(idelete)]
        umat=umat[:,:-len(idelete)]
    
    #Now go one epoch at a time, and interpolate all fields
    #Only apply topo effect to some points
    if coast_file!=None:
        #Straight line coordinates
        coast=genfromtxt(coast_file)
        #Get mask for applying horizontal effect
        mask=zeros(loni.shape)
        for k1 in range(loni.shape[0]):
            for k2 in range(loni.shape[1]):
                ip=argmin(abs(lati[k1,k2]-coast[:,1]))
                if loni[k1,k2]>coast[ip,0]:#Change this depending on which direction you want effect to be applied
                    print('Applying coastal mask')
                    mask[k1,k2]=nan
        imask1,imask2=where(mask==0)#Points tot he right DO apply horiz. effect
    print('... interpolating coseismic offsets to a regular grid')
    nt_iter=umat.shape[0]
    for kt in range(nt_iter):
        if kt%20==0:
            print('... ... working on time slice '+str(kt)+' of '+str(nt_iter))
        ninterp=griddata((lon,lat),nmat[kt,:],(loni,lati),method='cubic',fill_value=0)
        einterp=griddata((lon,lat),emat[kt,:],(loni,lati),method='cubic',fill_value=0)
        uinterp=griddata((lon,lat),umat[kt,:],(loni,lati),method='cubic',fill_value=0)
        #Output vertical
        uout=uinterp.copy()
        #Apply effect of topography advection
        if topo_effect==False:
            print('WARNING: No topography effect added')
        else:
            print('Applying topo effect')
            if coast_file==None: #Apply everywhere
                uout=uout+zdx*einterp+zdy*ninterp
            else: #Don't apply topo effect on dry land
                uout[imask1,imask2]=uout[imask1,imask2]+zdx[imask1,imask2]*einterp[imask1,imask2]+zdy[imask1,imask2]*ninterp[imask1,imask2]
        #print 'no horiz'
        #Filter?
        if variance!=None:
            uout=gaussian_filter(uout,variance)
        #Convert to column format and append
        xyz=grd2xyz(uout,loni,lati)
        tvec=(kt*tsun_dt)*ones((len(xyz),1))
        if kt==0: #Intialize
            numel=uout.size #Number of elements in grid
            kwrite=numel #Where to write the data
            dtopo=zeros((numel*nt_iter,4))
            dtopo[0:kwrite,1:3]=xyz[:,0:2]
        else:
            dtopo[kwrite:kwrite+numel,:]=c_[tvec,xyz]
            kwrite=kwrite+numel
        if static==True:
            tvec=ones(tvec.shape)
            numel=uout.size*2 #Number of elements in grid
            kwrite=numel/2 #Where to write the data
            dtopo=zeros((numel*nt_iter,4))
            dtopo[0:kwrite,1:3]=xyz[:,0:2]
            dtopo[kwrite:kwrite+numel,1:3]=xyz[:,0:2]
            dtopo[kwrite:kwrite+numel,:]=c_[tvec,xyz]
            kwrite=kwrite+numel
    print('... writting dtopo files')
    savetxt(data_dir+outname+'.dtopo',dtopo,fmt='%i\t%.6f\t%.6f\t%.4e')   
    print('Output to '+data_dir+outname+'.dtopo') 
            

def move_seafloor_okada(mudpy_file,out_file,x,y,refine_factor=None,mu=40e9,return_object=False):
    """
    Use the Okada routine in GeoClaw to generate the dtopo file
    
    x=arange(-75,-73,0.05)
    y=arange(-44.5,-42.5,0.05)
    """
    
    from clawpack.geoclaw import dtopotools
    from numpy import genfromtxt,arctan2,rad2deg,squeeze
    from copy import deepcopy
    
    #load rupture info
    f=genfromtxt(mudpy_file)
    
    #init object
    fault = dtopotools.Fault() 
              
    #Make one subfault object at a time
    for k in range(0,len(f)):
        subfault = dtopotools.SubFault()
        subfault.strike = f[k,4]
        subfault.dip = f[k,5]
        #Vertical faults introduce errors
        if subfault.dip==90:
            subfault.dip=89.9
        subfault.rake = rad2deg(arctan2(f[k,9],f[k,8]))
        subfault.mu=mu
        subfault.coordinate_specification = "centroid"                
        subfault.length = f[k,10]
        subfault.width = f[k,11]
        subfault.depth = f[k,3]*1000
        subfault.slip = (f[k,9]**2+f[k,8]**2)**0.5
        subfault.longitude = f[k,1]
        subfault.latitude = f[k,2]
        if refine_factor!=None:
            #How many times are we doing it?
            for ksplit in range(int(refine_factor)):
                if ksplit==0:    
                    subfault_list=split_subfault(subfault)
                else:
                    subfault_list=split_subfault(subfault_list)#len(split_subfault(subfault_list)) len(subfault_list)
        #Decide how to append     
        if k==0:
            if refine_factor==None:
                fault.subfaults=[subfault]
            else:
                fault.subfaults=deepcopy(subfault_list)
        else:
            if refine_factor==None:
                fault.subfaults.append(subfault)
            else:
                fault.subfaults.extend(subfault_list)
     
    ##Run Okada   
    fault.create_dtopography(x,y,times=[0.,1.],verbose=True)
    
    #Get dtopo
    dtopo=fault.dtopo
    dtopo.write(path=out_file,dtopo_type=1)
    
    if return_object:
        return fault
    


def split_subfault(subfault_list):
    '''
    Split a subfault object into 4 samller subfaults
    '''
    
    from pyproj import Geod
    from numpy import deg2rad,cos,sin
    from copy import deepcopy
    
    try:
        Nsubs=len(subfault_list)
    except:
        Nsubs=1
    #Loop over all elements in subfault_list
    for ksub in range(Nsubs):
        if Nsubs==1:
            subfault=deepcopy(subfault_list)
        else:
            subfault=deepcopy(subfault_list[ksub])
            
    
        #Assing useful object properties
        lon=subfault.longitude
        lat=subfault.latitude
        depth=subfault.depth
        strike=subfault.strike
        dip=subfault.dip
        L=subfault.length
        W=subfault.width
        
        #Projection object
        P=Geod(ellps='WGS84')
        
        #Lon lat of pseudo centers
        lon_P1,lat_P1,baz=P.fwd(lon,lat,strike,L/4.0)
        lon_P2,lat_P2,baz=P.fwd(lon,lat,strike,-L/4.0)
        
        # Pseudo widths
        wH=W*cos(deg2rad(dip))
        wV=W*sin(deg2rad(dip))
        
        #Get final lon lat depth of 4 new faults
        lon1,lat1,baz=P.fwd(lon_P1,lat_P1,strike-90,wH/4.0)
        lon2,lat2,baz=P.fwd(lon_P1,lat_P1,strike-90,-wH/4.0)
        z1=depth-wV/4.
        z2=depth+wV/4.
        
        lon3,lat3,baz=P.fwd(lon_P2,lat_P2,strike-90,wH/4.0)
        lon4,lat4,baz=P.fwd(lon_P2,lat_P2,strike-90,-wH/4.0)
        z3=depth-wV/4.
        z4=depth+wV/4.
        
        #Assemble into subfault objects and into list
        subfault1 = deepcopy(subfault)
        subfault2 = deepcopy(subfault)
        subfault3 = deepcopy(subfault)
        subfault4 = deepcopy(subfault)
        
        subfault1.length = subfault2.length = subfault3.length = subfault4.length = L/2.
        subfault1.width = subfault2.width = subfault3.width = subfault4.width = W/2.
        
        subfault1.longitude=lon1 ; subfault1.latitude=lat1 ; subfault1.depth=z1
        subfault2.longitude=lon2 ; subfault2.latitude=lat2 ; subfault2.depth=z2
        subfault3.longitude=lon3 ; subfault3.latitude=lat3 ; subfault3.depth=z3
        subfault4.longitude=lon4 ; subfault4.latitude=lat4 ; subfault4.depth=z4
        
        if ksub==0:
            subfault_list_out=[subfault1,subfault2,subfault3,subfault4]
        else:
            subfault_list_out.extend([subfault1,subfault2,subfault3,subfault4])
            
    return subfault_list_out




###########                Tools and trinkets                      #############
    
def get_mu(structure,zs,return_speeds=False):
    '''
    Look in structure model and compute rigidity given a source depth
    
    IN:
        structure: Array with velocity structure information
        zs: depth in km at which you want to compute mu
        
    OUT:
        mu: Rigidity in Pa
    '''
    from numpy import nonzero
    
    if len(structure)>1: #Model is more than jsut the halfspace
        Z=structure[:,0].cumsum()
        #Which layer do I want?
        i=nonzero(zs>Z)[0]
        if i.size==0: #It's in top layer
            imu=0
        else:
            imu=max(i)+1
        if imu>=structure.shape[0]:
            imu=imu-1#It's int he half-space
        mu=((1000*structure[imu,1])**2)*structure[imu,3]*1000
        beta=1000*structure[imu,1]
        alpha=1000*structure[imu,2]
        #print "Rigidity at z="+str(zs)+' is, mu = '+str(mu/1e9)+'GPa'
    else: #Model is a halfspace
        mu=((1000*structure[0,1])**2)*structure[0,3]*1000
        beta=1000*structure[0,1]
        alpha=1000*structure[0,2]
    
    if return_speeds==False:
        return mu
    else: 
        return mu,alpha,beta
        
        
def get_Q(structure,zs,perturb=0.0001):
    '''
    Look in structure model return Qp,Qs
    
    IN:
        structure: Array with velocity structure information
        zs: depth in km at which you want to compute mu
        
    OUT:
        Qp and Qs
    '''
    from numpy import nonzero
    
    if len(structure)>1: #Model is more than jsut the halfspace
        Z=structure[:,0].cumsum()
        #Which layer do I want?
        i=nonzero(zs>Z)[0]
        if i.size==0: #It's in top layer
            imu=0
        else:
            imu=max(i)+1
        if imu>=structure.shape[0]:
            imu=imu-1#It's int he half-space
        Qp=structure[imu,5]
        Qs=structure[imu,4]
    else: #Model is a halfspace
        Qp=structure[0,5]
        Qs=structure[0,4]
    
    return Qp,Qs
        
    
def get_mu_and_area(home,project_name,fault_name,model_name):
    
    from numpy import genfromtxt,zeros
    
    #load fault
    fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    #load structure
    structure=genfromtxt(home+project_name+'/structure/'+model_name)
    
    #get rigidities
    mu=zeros(len(fault))
    for k in range(len(mu)):
        depth=fault[k,3]
        mu[k]=get_mu(structure,depth)
        
    #Get areas
    area=fault[:,-2]*fault[:,-1]
    
    return mu,area
    

def get_source_time_function(mu,area,rise_time,t0,slip):
    '''
    Compute source time function for a given rise time, right now it assumes 1m of slip
    and a triangle STF
    '''
    from numpy import zeros,linspace,where
    
    rise_time=float(rise_time)
    #Initialize outputs
    t=linspace(t0,t0+rise_time,1000)
    Mdot=zeros(t.shape)
    #Triangle gradient
    m=4*mu*area/(rise_time**2)
    #Upwards intercept
    b1=-m*t0
    #Downwards intercept
    b2=m*(t0+rise_time)
    #Assign moment rate
    i=where(t<=t0+rise_time/2)[0]
    Mdot[i]=m*t[i]+b1
    i=where(t>t0+rise_time/2)[0]
    Mdot[i]=-m*t[i]+b2 
    Mdot=Mdot*slip  
    return t,Mdot
    
def add2stf(t1,Mdot1,t2,Mdot2):
    '''
    Add two overlapping source time functions
    '''
    from numpy import interp,linspace
    #Make interpolation vector
    tstart=min(t1[0],t2[0])
    tend=max(t1[-1],t2[-1])
    ti=linspace(tstart,tend,10000)
    #Interpolate
    Mdot1_interp=interp(ti,t1,Mdot1,left=0,right=0)
    Mdot2_interp=interp(ti,t2,Mdot2,left=0,right=0)
    #Add them up
    Mdot_out=Mdot1_interp+Mdot2_interp
    return ti,Mdot_out
    

def add_traces(ss,ds,ssmult,dsmult):
    '''
    Add two stream objects with dip slip and strike slip contributions. This code will take
    two stream objects and super impsoe them according tot he weights defined by ssmult
    dsmult. If one waveform is longer than the other then the code will extend the
    shorter waveform by padding it with the last value.
    
    For simple addition use ss=ds=M=1
    
    IN:
        ss: Strike slip waveform
        ds: Dip-slip waveform
        ssmult: Strike-slip contribution (meters)
        dsmult: Strike-slip contribution (meters)
    
    OUT:
        st: Stream object with result of superposition
    '''
    from numpy import zeros
    
    #If one stream object is empty set it to zeros, if both are empty then freak out
    if ss.count()==0 and ds.count()==0:
        print('FATAL ERROR: can\'t add 2 empty stream objects doofus')
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
    st[0].data=ss[0].data+ds[0].data
    #And done
    return st


def tshift(st,tshift):
    '''
    Shift a stream object by tshift seconds, positive moves forward in time
    
    IN:
        st: Stream object
        tshift: Number fo seconds to shift
    OUT:
        st: Shifted stream object
    '''
    from datetime import timedelta
    td=timedelta(seconds=tshift)
    st[0].stats.starttime=st[0].stats.starttime+td
    return st       
    
                        
def tshift_trace(nss,ess,zss,nds,eds,zds,tshift,time_epi,npts):
    '''
    Shift a stream object by tshift seconds, positive moves forward in time and then trim it
    so that it spans from time_epi to time_epi+(npts*dt).
    
    IN:
        st: Stream object
        tshift: Number fo seconds to shift
    OUT:
        st: Shifted stream object
    '''
    
    from datetime import timedelta
    from numpy import r_,zeros,ones
    
    # Shifted onset time
    tstart=nss.stats.starttime+timedelta(seconds=tshift)
    #How many positions from start of tr to start of trace?
    dt=nss.stats.delta
    nshift=int(round((tstart-time_epi)/dt))
    if nshift>0:
        #Place in output variable
        nss.data=r_[zeros(nshift),nss.data[0:npts-nshift]]
        ess.data=r_[zeros(nshift),ess.data[0:npts-nshift]]
        zss.data=r_[zeros(nshift),zss.data[0:npts-nshift]]
        nds.data=r_[zeros(nshift),nds.data[0:npts-nshift]]
        eds.data=r_[zeros(nshift),eds.data[0:npts-nshift]]
        zds.data=r_[zeros(nshift),zds.data[0:npts-nshift]]    
    else:
        nss.data=r_[nss.data[-nshift:],nss.data[-1]*ones(-nshift)]
        ess.data=r_[ess.data[-nshift:],ess.data[-1]*ones(-nshift)]
        zss.data=r_[zss.data[-nshift:],zss.data[-1]*ones(-nshift)]
        nds.data=r_[nds.data[-nshift:],nds.data[-1]*ones(-nshift)]
        eds.data=r_[eds.data[-nshift:],eds.data[-1]*ones(-nshift)]
        zds.data=r_[zds.data[-nshift:],zds.data[-1]*ones(-nshift)]
    
    return nss,ess,zss,nds,eds,zds
    
    


def trim_to_hypo_time_and_pad(tr,hypo_time,npts):
    """
    
    This function si tobe used for making sure syntehtics start at hypo time 
    and have the correct lengths. Will zero pad at begining or pad with last
    value at end of trace depedning on the situation

    Parameters
    ----------
    tr : obspy trace objet
        The inputt ace to be trimmed
    hypo_time : UTCDateTime object
        Hypocentral time for the rupture
    npts : int
        Total number of points expected in output trace, usually shoudl be 
        same as NFFT used in GFs/synthetics calculation

    Returns
    -------
    tr : obspy trace object
        Trimmed and paddedtrace

    """
    
    
    from numpy import zeros,r_,ones
    
    samples_to_pad = abs(int(round((tr.stats.starttime - hypo_time)/tr.stats.delta)))
    
    if samples_to_pad > npts:  #trace is waaaay too far into the future
        tr.data = zeros(npts)
        
    else:
        if tr.stats.starttime > hypo_time:
            tr.data = r_[zeros(samples_to_pad),tr.data[0:npts-samples_to_pad]]
    
        else:
            tr.data = r_[tr.data[samples_to_pad:],tr.data[-1]*ones(samples_to_pad)]
        
    #And update starttime
    tr.stats.starttime = hypo_time
    
    return tr    
    
    
                                                            
        
def round_time(t1,delta):
    '''
    Round the initial time of a waveform to the nearest multiple of the sampling rate
    IN:
        t1: UTC time object containing start time of waveform
        delta: Sampling interval of waveform in seconds
    OUT:
        t1: Rounded UTC time object
    '''
    from datetime import timedelta
    #Move start and end times to start exactly on the dt intervals
    dtmicro=delta*1e6
    intervals=t1.microsecond/dtmicro #How many intervals in microseconds
    adjustment=(round(intervals)-intervals)*dtmicro
    td=timedelta(microseconds=adjustment)
    t1=t1+td
    return t1

def upsample(st,delta):
    '''
    Go from a low sampling rate to a high sampling rate
    
    IN:
        st - stream object
        delta - sampling rate requested in seconds
    
    OUT:
        st - modified stream object
    '''
    
    from scipy.interpolate import interp1d
    from numpy import arange
    
    t=st[0].times()
    y=st[0].data
    #Make interpolant
    f=interp1d(t,y)
    ti=arange(t[0],t[-1],delta)
    #Interpoalte and reassign tos tream object
    yi=f(ti)
    st[0].data=yi
    st[0].stats.delta=delta

def lowpass(data,fcorner,fsample,order,zerophase=True,cascading=True):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt,lfilter,sosfiltfilt,sosfilt
    from numpy import size,array
    
    if size(fcorner)==2:
        ftype='bandpass'
    else:
        ftype='lowpass'
    fnyquist=fsample/2
    
    if cascading==True:
        sos = butter(order, array(fcorner)/(fnyquist),ftype,output='sos')
        if zerophase==True:
            data_filt=sosfiltfilt(sos,data)
        else:
            data_filt=sosfilt(sos,data)
    
    else:
        b, a = butter(order, array(fcorner)/(fnyquist),ftype)
        if zerophase==True:
            data_filt=filtfilt(b,a,data)
        else:
            data_filt=lfilter(b,a,data)
        
    return data_filt
    
def highpass(data,fcorner,fsample,order,zerophase=True,cascading=True):
    '''
    Make a highpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt,lfilter,sosfiltfilt,sosfilt
    from numpy import size,array
    
    fnyquist=fsample/2
    
    if cascading==True:
        sos = butter(order, array(fcorner)/(fnyquist),'highpass',output='sos')
        if zerophase==True:
            data_filt=sosfiltfilt(sos,data)
        else:
            data_filt=sosfilt(sos,data)
    
    else:
        b, a = butter(order, array(fcorner)/(fnyquist),'highpass')
        if zerophase==True:
            data_filt=filtfilt(b,a,data)
        else:
            data_filt=lfilter(b,a,data)
    
    return data_filt
 
def bandstop(data,fcorner,fsample,order,zerophase=True,cascading=True):
    '''
    Make a lowpass zero phase filter
    '''
    from scipy.signal import butter,filtfilt,lfilter,sosfiltfilt,sosfilt
    from numpy import array
    

    fnyquist=fsample/2
    
    if cascading==True:
        sos = butter(order, array(fcorner)/(fnyquist),'bandstop',output='sos')
        if zerophase==True:
            data_filt=sosfiltfilt(sos,data)
        else:
            data_filt=sosfilt(sos,data)
    
    else:
        b, a = butter(order, array(fcorner)/(fnyquist),'bandstop')
        if zerophase==True:
            data_filt=filtfilt(b,a,data)
        else:
            data_filt=lfilter(b,a,data)
        
    return data_filt   

    
def inv2coulomb(rupt,epicenter,fout,zeroslip=False,offset=0):
    '''
    Convert .inv file to Coulomb-ready .inp file
    
    IN:
        rupt - path ro rupture (.inv) file
    '''
    import pyproj
    from numpy import genfromtxt,unique,zeros,where,deg2rad,sin,cos
    
    #Read fault
    f=genfromtxt(rupt)
    #Get total slip by identifying unique fault numbers
    if zeroslip==False: # Get actaul slip
        u=unique(f[:,0])
        ss=zeros(len(u))
        ds=zeros(len(u))
        all_ss=f[:,8]
        all_ds=f[:,9]
        for k in range(len(u)):
            i=where(u[k]==f[:,0])
            ss[k]=all_ss[i].sum()
            ds[k]=all_ds[i].sum()
        #Sum them
        slip=(ss**2+ds**2)**0.5
        #Get rake
        rake=ssds2rake(ss,ds)
        Nfaults=len(u)
        #Get width and length to get coordiantes of top corners and strike,dip
        width=f[:,10]/1000
        length=f[:,11]/1000
    else: #Make zero slip receiver faults
        slip=zeros(len(f))
        ss=zeros(len(f))
        ds=zeros(len(f))
        Nfaults=len(f)
        width=f[:,8]/1000
        length=f[:,9]/1000
    #Convert coordinate subfault centers to local cartesian
    g = pyproj.Geod(ellps='WGS84') # Use WGS84 ellipsoid.
    x=zeros(Nfaults)
    y=zeros(Nfaults)
    for k in range(Nfaults):
        baz,az,d=pyproj.Geod.inv(g,f[k,1],f[k,2],epicenter[0],epicenter[1])
        x[k]=(d/1000)*sin(deg2rad(az))
        y[k]=(d/1000)*cos(deg2rad(az))

    strike=f[:,4]
    dip=f[:,5]
    depth=f[:,3]
    top_mid_x=zeros(Nfaults)
    top_mid_y=zeros(Nfaults)
    top_direction=strike-90 #This is the angle that points towards the top edge of the fault
    xstart=zeros(Nfaults)
    ystart=zeros(Nfaults)
    xfin=zeros(Nfaults)
    yfin=zeros(Nfaults)
    ztop=zeros(Nfaults)
    zbot=zeros(Nfaults)
    for k in range(Nfaults):
        top_mid_x[k]=x[k]+((width[k]/2)*cos(deg2rad(dip[k])))*sin(deg2rad(top_direction[k]))
        top_mid_y[k]=y[k]+((width[k]/2)*cos(deg2rad(dip[k])))*cos(deg2rad(top_direction[k]))
        xstart[k]=top_mid_x[k]+(width[k]/2)*sin(deg2rad(strike[k]-180))
        ystart[k]=top_mid_y[k]+(width[k]/2)*cos(deg2rad(strike[k]-180))
        xfin[k]=top_mid_x[k]+(width[k]/2)*sin(deg2rad(strike[k]))
        yfin[k]=top_mid_y[k]+(width[k]/2)*cos(deg2rad(strike[k]))
        ztop[k]=depth[k]-(length[k]/2)*sin(deg2rad(dip[k]))
        zbot[k]=depth[k]+(length[k]/2)*sin(deg2rad(dip[k]))
    #Write to file and then manually add the ehaders and footers by copy pasting from some NEIC file (LAZY!)
    f=open(fout,'w')
    for k in range(Nfaults):
        #out='1   %10.4f %10.4f %10.4f %10.4f 100 %10.4f %10.4f %10.4f %10.4f %10.4f %i\n' % (xstart[k],ystart[k],xfin[k],yfin[k],rake[k],slip[k],dip[k],ztop[k],zbot[k],k)
        out='%3d %10.4f %10.4f %10.4f %10.4f 100 %10.4f %10.4f %10.4f %10.4f %10.4f %d\n' % (1,xstart[k],ystart[k],xfin[k],yfin[k],-ss[k],ds[k],dip[k],ztop[k],zbot[k],k+1+offset)
        f.write(out)
    f.close()
    
def coulomb_xy2latlon(f,epicenter,fout):
    '''
    Change the x-y coordinates of a Coulomb file to lat/lon
    '''
    from numpy import genfromtxt,zeros,rad2deg,arctan,isnan,savetxt
    import pyproj
    
    s=genfromtxt(f)
    x=s[:,1]
    y=s[:,2]
    #Now use pyproj to dead reckon anf get lat/lon coordinates of subfaults
    g = pyproj.Geod(ellps='WGS84')
    #first get azimuths of all points, go by quadrant
    az=zeros(x.shape)
    for k in range(len(x)):
        if x[k]>0 and y[k]>0:
            az[k]=rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]>0:
            az[k]=360+rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
        if x[k]>0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
    #Quadrant correction
    #Now horizontal distances
    d=((x**2+y**2)**0.5)*1000
    #Now reckon
    lo=zeros(len(d))
    la=zeros(len(d))
    for k in range(len(d)):
        if isnan(az[k]): #No azimuth because I'm on the epicenter
            print('Point on epicenter')
            lo[k]=epicenter[0]
            la[k]=epicenter[1]
        else:
            lo[k],la[k],ba=g.fwd(epicenter[0],epicenter[1],az[k],d[k])
    s[:,1]=lo
    s[:,2]=la
    savetxt(fout,s,fmt='%.6f')
    
def coulomb_disp_xy2latlon(f,epicenter,fout):
    '''
    Change the x-y coordinates of a Coulomb file to lat/lon
    '''
    from numpy import genfromtxt,zeros,rad2deg,arctan,isnan,savetxt,where
    import pyproj
    
    s=genfromtxt(f)
    x=s[:,0]
    y=s[:,1]
    #Fix problem with x=0 and y=0 lines
    i=where(x==0)[0]
    x[i]=0.001
    i=where(y==0)[0]
    y[i]=0.001
    #Now use pyproj to dead reckon anf get lat/lon coordinates of subfaults
    g = pyproj.Geod(ellps='WGS84')
    #first get azimuths of all points, go by quadrant
    az=zeros(x.shape)
    for k in range(len(x)):
        if x[k]>0 and y[k]>0:
            az[k]=rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]>0:
            az[k]=360+rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
        if x[k]>0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
    #Quadrant correction
    #Now horizontal distances
    d=((x**2+y**2)**0.5)*1000
    #Now reckon
    lo=zeros(len(d))
    la=zeros(len(d))
    for k in range(len(d)):
        if isnan(az[k]): #No azimuth because I'm on the epicenter
            print('Point on epicenter')
            lo[k]=epicenter[0]
            la[k]=epicenter[1]
        else:
            lo[k],la[k],ba=g.fwd(epicenter[0],epicenter[1],az[k],d[k])
    s[:,0]=lo
    s[:,1]=la
    savetxt(fout,s,fmt='%.6f')
    
      
def ssds2rake(ss,ds):
    '''
    Compute rake angle in degrees
    '''
    from numpy import arctan,rad2deg,pi,zeros
    try:
        rake=zeros(len(ss))
    except:
        rake=zeros(1)
    for k in range(len(rake)):
        if ss[k]>0 and ds[k]>0:
            rake[k]=arctan(ds[k]/ss[k])
        elif ds[k]>0 and ss[k]<0:
            rake[k]=pi+arctan(ds[k]/ss[k])
        elif ds[k]<0 and ss[k]<0:
            rake[k]=pi+arctan(ds[k]/ss[k])
        elif ds[k]<0 and ss[k]>0:
            rake[k]=2*pi+arctan(ds[k]/ss[k])
    rake=rad2deg(rake)
    return rake
    
def makefault(fout,strike,dip,nstrike,dx_dip,dx_strike,epicenter,num_updip,num_downdip,rise_time):
    '''
    Make a planar fault
    
    strike - Strike angle (degs)
    dip - Dip angle (degs)200/5
    '''
    from numpy import arange,sin,cos,deg2rad,r_,ones,arctan,rad2deg,zeros,isnan,unique,where,argsort
    import pyproj
    
    proj_angle=180-strike #Angle to use for sin.cos projection (comes from strike)
    y=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    x=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    z=ones(x.shape)*epicenter[2]
    y=y*cos(deg2rad(strike))
    x=x*sin(deg2rad(strike))   
    #Save teh zero line
    y0=y.copy()
    x0=x.copy()
    z0=z.copy()
    #Initlaize temp for projection up/down dip
    xtemp=x0.copy()
    ytemp=y0.copy()
    ztemp=z0.copy()
    #Get delta h and delta z for up/ddx_dip=1own dip projection
    if num_downdip>0 and num_updip>0:
        dh=dx_dip*cos(deg2rad(dip))
        dz=dx_dip*sin(deg2rad(dip))
        #Project updip lines
        for k in range(num_updip):
            xtemp=xtemp+dh*cos(deg2rad(proj_angle))
            ytemp=ytemp+dh*sin(deg2rad(proj_angle))
            ztemp=ztemp-dz
            x=r_[x,xtemp]
            y=r_[y,ytemp]
            z=r_[z,ztemp]
        #Now downdip lines
        xtemp=x0.copy()
        ytemp=y0.copy()
        ztemp=z0.copy()
        for k in range(num_downdip):
            xtemp=xtemp-dh*cos(deg2rad(proj_angle))
            ytemp=ytemp-dh*sin(deg2rad(proj_angle))
            ztemp=ztemp+dz
            x=r_[x,xtemp]
            y=r_[y,ytemp]
            z=r_[z,ztemp]
    #Now use pyproj to dead reckon anf get lat/lon coordinates of subfaults
    g = pyproj.Geod(ellps='WGS84')
    #first get azimuths of all points, go by quadrant
    az=zeros(x.shape)
    for k in range(len(x)):
        if x[k]>0 and y[k]>0:
            az[k]=rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]>0:
            az[k]=360+rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
        if x[k]>0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
    #Quadrant correction
    #Now horizontal distances
    d=((x**2+y**2)**0.5)*1000
    #Now reckon
    lo=zeros(len(d))
    la=zeros(len(d))
    for k in range(len(d)):
        if isnan(az[k]): #No azimuth because I'm on the epicenter
            print('Point on epicenter')
            lo[k]=epicenter[0]
            la[k]=epicenter[1]
        else:
            lo[k],la[k],ba=g.fwd(epicenter[0],epicenter[1],az[k],d[k]) 
    #Sort them from top right to left along dip
    zunique=unique(z)
    for k in range(len(zunique)):
        i=where(z==zunique[k])[0] #This finds all faults at a certain depth
        isort=argsort(la[i]) #This sorths them south to north
        if k==0: #First loop
            laout=la[i][isort]
            loout=lo[i][isort]
            zout=z[i][isort]
        else:
            laout=r_[laout,la[i][isort]]
            loout=r_[loout,lo[i][isort]]
            zout=r_[zout,z[i][isort]]
    #Write to file
    strike=ones(loout.shape)*strike
    dip=ones(loout.shape)*dip
    tw=ones(loout.shape)*0.5
    rise=ones(loout.shape)*rise_time
    L=ones(loout.shape)*dx_strike*1000
    W=ones(loout.shape)*dx_dip*1000
    f=open(fout,'w')
    for k in range(len(x)):   
        out='%i\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\n' % (k+1,loout[k],laout[k],zout[k],strike[k],dip[k],tw[k],rise[k],L[k],W[k])
        f.write(out)
    f.close()
    
    
def make_checkerboard(rupt,nstrike,ndip,fout,nwin):
    '''
    Make checekrboard for resolution analysis
    '''
    from numpy import genfromtxt,unique,arange,union1d,savetxt
   
    f=genfromtxt(rupt)
    u=unique(f[:,0])
    f=f[0:len(u),:]
    #Set strike-slip/dip_slip to zero
    f[:,8]=0
    f[:,9]=0
    #1 and 2
    i1=arange(0,nstrike,4)
    i2=arange(1,nstrike,4)
    i=union1d(i1,i2)
    j1=arange(21,2*nstrike,4)
    j2=arange(22,2*nstrike,4)
    j=union1d(j1,j2)
    ij=union1d(i,j)
    f[ij,9]=1
    # 3 and 4
    i1=arange(44,3*nstrike,4)
    i2=arange(45,3*nstrike,4)
    i=union1d(i1,i2)
    j1=arange(65,4*nstrike,4)
    j2=arange(66,4*nstrike,4)
    j=union1d(j1,j2)
    ij=union1d(i,j)
    f[ij,9]=1
    # 5 and 6
    i1=arange(84,5*nstrike,4)
    i2=arange(85,5*nstrike,4)
    i=union1d(i1,i2)
    j1=arange(105,6*nstrike,4)
    j2=arange(106,6*nstrike,4)
    j=union1d(j1,j2)
    ij=union1d(i,j)
    f[ij,9]=1
    # 7 and 8
    i1=arange(128,7*nstrike,4)
    i2=arange(129,7*nstrike,4)
    i=union1d(i1,i2)
    j1=arange(149,8*nstrike,4)
    j2=arange(150,8*nstrike,4)
    j=union1d(j1,j2)
    ij=union1d(i,j)
    f[ij,9]=1
    # 9
    i1=arange(168,9*nstrike,4)
    i2=arange(169,9*nstrike,4)
    ij=union1d(i1,i2)
    f[ij,9]=1
    #Write to file
    fmt='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
    savetxt(fout,f,fmt)
    
def model_resolution(home,project_name,GF_list,Gfile,Ltype,lambda_s,nstrike,ndip,bandpass,G=None,Ls=None,G_from_file=False,return_Ls=False):
    '''
    Compute model resolution matrix and output to GMT plottable file
    '''
    from numpy import load,diag
    from scipy.linalg import inv,norm
    from mudpy.inverse import getdata,get_data_weights
    
    #load data vector (needed for normalization weight)
    decimate=None
    d=getdata(home,project_name,GF_list,decimate,bandpass=bandpass)
    
    #get weights used in inversion
    w=get_data_weights(home,project_name,GF_list,d,decimate)
    
    #apply to data vector
    d=w*d.squeeze()
    
    #normalize
    data_norm = norm(d)
    d = d/data_norm
    
    #Load Green's functions
    if G_from_file: #load from file
        G=load(Gfile)
    else: #Use G from user
        pass
    
    #weight G according tow eights used in ivnersion
    W=diag(w)
    
    #Normalization effect
    W = W/data_norm
    G=W.dot(G)
    
    if Ls is None: #Make regularization matrix
        top='free' ; bottom='locked' ; left='locked' ; right='locked'
        bounds=(top,bottom,left,right)
        Ls = laplacian_smoother(nstrike,ndip,bounds)
    else: #Use provided by user
        pass
    
    #Now compute the resolution matrix
    Gg = (G.T).dot(G) + (lambda_s**2)*(Ls.T).dot(Ls)
    Gg = inv(Gg).dot(G.T)
    
    R = Gg.dot(G) 
    
    if return_Ls == False:
        return R
    else:
        return R,Ls



def model_spread(R,fault,dist=None,return_dist=False):
    '''
    Compute model resolution spread
    dist is a matrix of itnerfault distances
    '''

    from pyproj import Geod
    from numpy import zeros
    
    #projection object
    p=Geod(ellps='WGS84')
    
    if dist is None: #calcualte dist matrix
        #build distances matrix
        dist=zeros(R.shape)
        for row in range(R.shape[0]):
            for col in range(R.shape[1]):
                
                sub1=fault[row,1:4]
                sub2=fault[col,1:4]
                az,baz,dist_in_km=p.inv(sub1[0],sub1[1],sub2[0],sub2[1])
                dist_in_km /= 1000 # to km
                
                #add depth
                dist_in_km = (dist_in_km**2+(sub1[2]-sub2[2])**2)**0.5
                dist[row,col] = dist_in_km
    else: #Just use the dist matrix provided by user
        pass            

    #calculate spread
    Rspread=dist*(R**2)
    spread=Rspread.sum()
    row_spread=Rspread.sum(axis=1)
    
    if return_dist==False:
        return spread,row_spread
    else:
        return spread,row_spread,dist
    
    
def laplacian_smoother(nstrike,ndip,bounds):
    '''
    Make spatial regularization matrix based on finite difference Lapalce operator.
    This routine will request adjustments depending on the boundary conditions requested
    on the edges of the fault model.
    '''
    
    from numpy import zeros
    from mudpy.inverse import laplace_stencil
    
    N=nstrike*ndip
    #Initalize
    L=zeros((2*N,2*N))
    #Which L am I building?
    print('Making discrete Laplace operator regularization matrix...')
    for kfault in range(N):#Loop over faults and fill regularization matrix
        print(kfault)
        print(nstrike)
        print(ndip)
        print(bounds)
        stencil,values=laplace_stencil(kfault,nstrike,ndip,bounds)
        #Add strike slip branches of stencil
        print('kfault:'+str(kfault))
        print('stencil:'+str(stencil))
        L[2*kfault,2*stencil]=values
        #Add dip slip branches of stencil
        L[2*kfault+1,2*stencil+1]=values
        Lout=L 
    return Lout
    
def trim_add_noise(data_path,checker_path,search_pattern):
    '''
    Trim checkerboard data and Add gaussian noise to data
    
    data_path='/Volumes/Kanagawa/Slip_Inv/tohoku_10s/data/waveforms/'
    search_pattern='checker.*disp*'
    checker_path='/Volumes/Kanagawa/Slip_Inv/tohoku_10s/output/forward_models/'    
    '''
    from numpy import var
    from numpy.random import normal
    from glob import glob
    from obspy import read
    
    checker_files=glob(checker_path+search_pattern)
    for k in range(len(checker_files)):
        ch=read(checker_files[k])
        #Find corresponding data file
        sta=checker_files[k].split('/')[-1].split('.')[1]
        vord=checker_files[k].split('/')[-1].split('.')[2]
        comp=checker_files[k].split('/')[-1].split('.')[3]
        data_file=glob(data_path+sta+'*'+vord+'*'+comp)
        st=read(data_file[0])
        ch.trim(starttime=st[0].stats.starttime,endtime=st[0].stats.endtime)
        #determine variance
        v=2e-5 #vel
        noise=normal(loc=0.0, scale=v**0.5, size=ch[0].stats.npts)
        ch[0].data=ch[0].data+noise
        ch.write(checker_files[k],format='SAC')
        
def padGFs(pad):
    '''
    Pad GFs with some extra zeros or with last value
    '''
    from glob import glob
    from obspy import read
    from numpy import ones,r_,zeros
    pad=150
    folders=glob('/Users/dmelgar/Slip_inv/Nepal_fwd/GFs/dynamic/*sub*')
    for k in range(len(folders)):
        print(str(k)+' / '+str(len(folders)))
        esubs=glob(folders[k]+'/*vel.e')
        nsubs=glob(folders[k]+'/*vel.n')
        zsubs=glob(folders[k]+'/*vel.z')
        for ksub in range(len(esubs)):
            e=read(esubs[ksub])
            n=read(nsubs[ksub])
            z=read(zsubs[ksub])
            #If displacement
            #e[0].data=r_[e[0].data,ones(pad)*e[0].data[-1]]
            #e.write(esubs[ksub],format='SAC')
            #
            #n[0].data=r_[n[0].data,ones(pad)*n[0].data[-1]]
            #n.write(nsubs[ksub],format='SAC')
            #
            #z[0].data=r_[z[0].data,ones(pad)*z[0].data[-1]]
            #z.write(zsubs[ksub],format='SAC')
            
            #If velocity
            e[0].data=r_[e[0].data,zeros(pad)]
            e.write(esubs[ksub],format='SAC')
            
            n[0].data=r_[n[0].data,zeros(pad)]
            n.write(nsubs[ksub],format='SAC')
            
            z[0].data=r_[z[0].data,zeros(pad)]
            z.write(zsubs[ksub],format='SAC')
            
def make_grid(lon_min,lon_max,lat_min,lat_max,delta_lon,delta_lat,out_file):
    '''
    Make a regular grid of points for GF computations
    '''
    
    from numpy import arange
    
    lon=arange(lon_min,lon_max+delta_lon,delta_lon)
    lat=arange(lat_min,lat_max+delta_lat,delta_lat)
    k=0
    f=open(out_file,'w')
    for i in range(len(lon)):
        for j in range(len(lat)):
            out='%s\t%10.4f\t%10.4f\n' %('SF'+str(k).rjust(4,'0'),lon[i],lat[j])
            f.write(out)
            k+=1
    f.close()
        
    
def usgs2fault(usgs_model,out_file,Dx,Dy,rise_time):
    '''
    Convert USGS finite fault to .fault
    '''
    from numpy import genfromtxt,ones,arange,savetxt,c_
    
    lon=genfromtxt(usgs_model,usecols=1)
    lat=genfromtxt(usgs_model,usecols=0)
    z=genfromtxt(usgs_model,usecols=2)
    st=genfromtxt(usgs_model,usecols=5)
    dip=genfromtxt(usgs_model,usecols=6)
    
    no=arange(1,len(lon)+1)
    H=Dx*ones(len(lon))
    W=Dy*ones(len(lon))
    tri=0.5*ones(len(lon))
    rt=rise_time*ones(len(lon))
    
    out=c_[no,lon,lat,z,st,dip,tri,rt,H,W]
    savetxt(out_file,out,fmt='%d\t%10.4f\t%10.4f\t%8.4f\t%6.1f\t%6.1f\t%.1f\t%.1f\t%.1f\t%.1f')
    
    
def usgs2rupt(usgs_model,out_file,Dx,Dy,N_header_lines= 9,rise_time_multiplier=1):
    '''
    Convert USGS finite fault to .fault
    '''
    from numpy import genfromtxt,ones,arange,savetxt,c_,deg2rad,cos,sin,zeros
    
    lon=genfromtxt(usgs_model,usecols=1,skip_header = N_header_lines)
    lat=genfromtxt(usgs_model,usecols=0,skip_header = N_header_lines)
    z=genfromtxt(usgs_model,usecols=2,skip_header = N_header_lines)
    st=genfromtxt(usgs_model,usecols=5,skip_header = N_header_lines)
    dip=genfromtxt(usgs_model,usecols=6,skip_header = N_header_lines)
    rake=genfromtxt(usgs_model,usecols=4,skip_header = N_header_lines)
    slip=genfromtxt(usgs_model,usecols=3,skip_header = N_header_lines)
    time=genfromtxt(usgs_model,usecols=7,skip_header = N_header_lines)
    
    tup = genfromtxt(usgs_model,usecols=8,skip_header = N_header_lines)
    tdown = genfromtxt(usgs_model,usecols=9,skip_header = N_header_lines)
    rise = tup + tdown
    fraction = tup / rise
    
    rise = rise * rise_time_multiplier
    #Parse rake out into SS and DS
    ss=(slip/100)*cos(deg2rad(rake))
    ds=(slip/100)*sin(deg2rad(rake))
    
    no=arange(1,len(lon)+1)
    H=Dx*ones(len(lon))
    W=Dy*ones(len(lon))
    mu=zeros(len(lon))

    out=c_[no,lon,lat,z,st,dip,fraction,rise,ss,ds,H,W,time,mu]
    savetxt(out_file,out,fmt='%d\t%10.4f\t%10.4f\t%8.4f\t%6.1f\t%6.1f\t%6.4f\t%6.2f\t%8.4f\t%8.4f\t%.1f\t%.1f\t%.2f\t%d')
    
def grid2xyz(home,project_name,sta_file,out_file):
    '''
    Read forward calculated station files and make an xyz file for plotting
    '''
    
    from numpy import genfromtxt
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+sta_file,usecols=0,dtype='U')
    lon=genfromtxt(home+project_name+'/data/station_info/'+sta_file,usecols=1)
    lat=genfromtxt(home+project_name+'/data/station_info/'+sta_file,usecols=2)
    f=open(out_file,'w')
    for k in range(len(sta)):
        neu=genfromtxt(home+project_name+'/output/forward_models/'+sta[k]+'.static.neu')
        n=neu[0]
        e=neu[1]
        u=neu[2]
        out='%s\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n' %(sta[k],lon[k],lat[k],n,e,u)
        f.write(out)
    f.close()
        
    
def static2kineamtic(rupt,epicenter,vr,fout):
    '''
    Take a static .inv or .rupt file and convert to a kinematic file given
    a rupture velocity vr
    '''
    from numpy import genfromtxt,savetxt
    from mudpy.inverse import epi2subfault
    
    #Read rupture file
    source=genfromtxt(rupt)
    #Determien distance from subfaults to epicenter
    tdelay=epi2subfault(epicenter,source,vr,0)
    source[:,12]=tdelay
    fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
    savetxt(fout,source,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')



def mudpy2sw4source(rupt,time_offset=0.0):
    '''
    Convert a mudpy rupt file to source type statememnts used by SW4 .in files
    '''
    
    
    from numpy import genfromtxt,sqrt,arctan2,rad2deg,array,pi
    from string import replace
    from os import path
    from obspy.imaging.scripts import mopad

    
    #Read mudpy file
    f=genfromtxt(rupt)

    #Define paths to output files
    folder=path.dirname(rupt)
    basename=path.basename(rupt)
    basename=replace(basename,'rupt','sw4source')
    fout=open(folder+'/'+basename,'w')
    
    #loop over subfaults
    for kfault in range(len(f)):
        
        zero_slip=False
        
        #Get subfault parameters
        lon=f[kfault,1]
        lat=f[kfault,2]
        depth=f[kfault,3]*1000 #in m for sw4
        strike=f[kfault,4]
        dip=f[kfault,5]
        area=f[kfault,10]*f[kfault,11] #in meters, cause this isn't dumb SRF
        tinit=f[kfault,12]+time_offset
        rake=rad2deg(arctan2(f[kfault,9],f[kfault,8]))
        slip=sqrt(f[kfault,8]**2+f[kfault,9]**2)
        rise_time=f[kfault,7]
        rigidity=f[kfault,13]

        #If subfault has zero rise time or zero slip
        zero_slip=False
        if slip==0:
            zero_slip=True
            print('Zero slip at '+str(kfault))
        elif rise_time==0:
            slip=0
            zero_slip=True
            print('Zero rise time at '+str(kfault))     
            
        #make rake be -180 to 180
        if rake>180:
            rake=rake-360
            
        #Make moment tensor for subfault (this is unit MT)
        MT=mopad.MomentTensor([strike,dip,rake])
        MT=array(MT.get_M())
        
        #Multiply by moment, remember M0=norm(MT)/sqrt(2)
        moment=rigidity*slip*area
        MT=MT*moment
        
        #Only write it if the thing has non-zero moment
        if zero_slip==False:
            mxx=MT[0,0]
            mxy=MT[0,1]
            mxz=MT[0,2]
            myy=MT[1,1]
            myz=MT[1,2]
            mzz=MT[2,2]
            freq=(2*pi)/rise_time
            # source lat=38.1035 lon=-120.8610 z=10000 mxx=-.0677e21 mxy=0.3149e21 mxz=0.2529e21 myy=-0.7636e21 myz=-0.5946e21 mzz=0.8313e21 type=Liu t0=3 freq=1
            line='source lat=%.6f lon=%.6f z=%.1f mxx=%.5e mxy=%.5e mxz=%.5e myy=%.5e myz=%.5e mzz=%.5e type=Liu t0=%.3f freq=%.3f\n' % (lat,lon,depth,mxx,mxy,mxz,myy,myz,mzz,tinit,freq)
            fout.write(line)
            
    fout.close()
            




def mudpy2srf(rupt,log_file,stf_dt=0.1,stf_type='triangle',inv_or_rupt='rupt',Ndip=None,
              hypocenter=None,time_pad=5.0,minSTFpoints=16,integrate=False,Mw=None):
    '''
    Convert a mudpy .rupt or .inv file to SRF version 2 format
    
    See for SRF format description here:
        https://scec.usc.edu/scecpedia/Standard_Rupture_Format
        
    if Ndip != None it means try to apportion Nstrike and Ndip to avoid huge Nstrike number
    '''
    
    from numpy import genfromtxt,where,sin,deg2rad,array,argmin,sqrt,sign,arctan2,rad2deg,arange,ones,zeros,r_
    from pyproj import Geod
    from os import path
    from scipy.integrate import cumtrapz
    
    #Define paths to output files
    folder=path.dirname(rupt)
    basename=path.basename(rupt)
    
    if inv_or_rupt == 'rupt':
        basename=basename.replace('rupt','srf')
        srf_file=folder+'/'+basename
        basename=basename.replace('srf','src')
        src_file=folder+'/'+basename
    else:
        basename=basename.replace('inv','srf')
        srf_file=folder+'/'+basename
        basename=basename.replace('srf','src')
        src_file=folder+'/'+basename       
    
    #Read mudpy file
    f=genfromtxt(rupt)
    
    
    #Calculate stuff for the header
    
    #Get coordiantes of shallowest row
    min_row=where(f[:,3]==f[:,3].min())[0]
    fmin=f[min_row,:]
    
    #Top center longitude
    elon=fmin[:,1].mean()
    
    #Top center latitude
    elat=fmin[:,2].mean()
    
    if Ndip==None:
        #Number of faults along strike
        Nstrike=len(fmin)
        #Number of faults along dip
        Ndip=len(f)/Nstrike
    else:
        #Number of faults along strike
        Nstrike=len(f)/Ndip
        
    #Subfault dimensions (in km)
    dx_strike=f[0,10]/1000.0
    dx_dip=f[0,11]/1000.0
    
    #Segment length and width
    length=Nstrike*dx_strike
    width=Ndip*dx_dip
    
    #Segment strike and dip
    strike=f[0,4]
    dip=f[0,5]
    
    #Depth to top of fault
    depth_top=fmin[0,3]-sin(deg2rad(dip))*(dx_dip/2)
    
    if inv_or_rupt == 'rupt':
        #Location of hypocenter from log file and magnitude
        flog=open(log_file,'r')
        while True:
            line=flog.readline()
            if 'Hypocenter (lon,lat,z[km])' in line:                
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
                hypocenter=array(s.split(',')).astype('float')
            if 'Actual magnitude' in line:
                Mw=float(line.split()[-1])
            elif line=='':
                break 
        flog.close()
    else:
        pass
    
    #Down dip location with regards to top of fault
    dip_hypo_position=(hypocenter[2]-depth_top)/sin(deg2rad(dip))
    
    #Get hypocenter row of subfaults and find coordinates of middle
    i=argmin(abs(f[:,3]-hypocenter[2]))
    test_depth=f[i,3]
    hypo_row=where(f[:,3]==test_depth)[0]
    hypo_center_lon=f[hypo_row,1].mean()
    hypo_center_lat=f[hypo_row,2].mean()
    #Distance from edge
    g=Geod(ellps='WGS84')
    az,baz,dist_from_center=g.inv(hypocenter[0],hypocenter[1],hypo_center_lon,hypo_center_lat)
    
    #Now decide if distance is positive or negative
    if strike>180:
        strike_rectified=strike-360
    else:
        strike_rectified=strike
    
    if sign(az)==sign(strike_rectified):
        strike_hypo_position=-dist_from_center/1000
    else:
        strike_hypo_position=dist_from_center/1000
        
    #open SRF file and write header data
    fout=open(srf_file,'w')
    fout.write('2.0\n') #SRF version
    fout.write('PLANE 1\n')
    fout.write('  %.6f\t%.6f\t%d\t%d\t%.4f\t%.4f\n' % (elon,elat,Nstrike,Ndip,length,width))
    fout.write('  %.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (strike,dip,depth_top,strike_hypo_position,dip_hypo_position))
    fout.write('POINTS %d\n' % (Nstrike*Ndip))
    
    #While we're here let's write the SRC file too
    fsrc=open(src_file,'w')
    fsrc.write('MAGNITUDE = %.4f\n' % Mw)
    fsrc.write('FAULT_LENGTH = %.4f\n' % length)
    fsrc.write('DLEN = %.4f\n' % dx_strike)
    fsrc.write('FAULT_WIDTH = %.4f\n' % width)
    fsrc.write('DWID = %.4f\n' % dx_dip)
    fsrc.write('DEPTH_TO_TOP = %.4f\n' % depth_top)
    fsrc.write('STRIKE = %.4f\n' % strike)
    fsrc.write('RAKE = 9999\n')
    fsrc.write('DIP = %.4f\n' % dip)
    fsrc.write('LAT_TOP_CENTER = %.4f\n' % elat)
    fsrc.write('LON_TOP_CENTER = %.4f\n' % elon)
    fsrc.write('HYPO_ALONG_STK = %.4f\n' % strike_hypo_position)
    fsrc.write('HYPO_DOWN_DIP = %.4f\n' % dip_hypo_position)
    fsrc.write('DT = 0.01\n')
    fsrc.write('SEED = 1')
    fsrc.close()
    
    #And that was jsut the headers, now let's move on to getting the subfault source time functions
    # Note mudpy works in mks SRF is cgs so must convert accordingly
    minNTstf=99999
    for kfault in range(len(f)):
        print(kfault)
        zero_slip=False
        #Get values for "Headers" for this subfault
        lon=f[kfault,1]
        lat=f[kfault,2]
        depth=f[kfault,3]
        strike=f[kfault,4]
        dip=f[kfault,5]
        area=f[kfault,10]*f[kfault,11]*100**2
        tinit=f[kfault,12]
        vs=2.80000e+05    #Default value for not known
        density=2.70000e+00 #default value for not known
        rake=rad2deg(arctan2(f[kfault,9],f[kfault,8]))
        slip=sqrt(f[kfault,8]**2+f[kfault,9]**2)*100
        
        #Now get source time function
        rise_time=f[kfault,7]
        total_time=rise_time*1.5
        #If subfault has zero rise time make it have tiny slip rate
        if slip==0:
            zero_slip=True
            stf = zeros(int(total_time/stf_dt))
            print('Zero slip at '+str(kfault))
        elif rise_time==0:
            slip=0
            zero_slip=True
            stf = zeros(int(total_time/stf_dt))
            print('Zero rise time at '+str(kfault))
        else:
            tstf,stf=build_source_time_function(rise_time,stf_dt,total_time,stf_type=stf_type,zeta=0.2,scale=True)
            #tstf,stf=build_source_time_function(rise_time,stf_dt,total_time,stf_type='triangle',scale=True)
            #Scale stf so integral matches total slip
            stf_adjust_factor=slip/stf_dt
            stf=stf*stf_adjust_factor #now tf is in cm/sec       
        
        #Now zero pad before and after end because SW4 struggles if subfault STFs are not zero padded
        if time_pad != None:
            zeros_pad=zeros(int(time_pad/stf_dt))
            stf=r_[zeros_pad,stf,zeros_pad]
            #Change start time of STF, it should now begin time_pad seconds earlier
            tinit=tinit-time_pad
        
        #How mant STF points?
        NTstf=len(stf)
        if NTstf<minSTFpoints: #Too short, zero pad
            print('Padding short STF...')
            zeros_pad=zeros(int(minSTFpoints/2))
            stf=r_[zeros_pad,stf,zeros_pad]
            #Change start time of STF, it should now begin time_pad seconds earlier
            time_shift=int(minSTFpoints/2)*stf_dt
            tinit=tinit-time_shift    
        
        
        #Check that everything is ok
        NTstf=len(stf)
        if NTstf<minNTstf:
            minNTstf=NTstf
            
        if zero_slip==True:
            NTstf=0
        
        #Write the subfault headers
        fout.write('  %.6f  %.6f  %.5e  %.2f  %.2f  %.5e  %.4f  %.4e  %.4e  %.4e\n' % (lon,lat,depth,strike,dip,area,tinit,stf_dt,vs,density))
        fout.write('  %.2f  %.4f  %d  0  0  0  0\n' % (rake,slip,NTstf))
            
        
        if zero_slip==False:
            #Integrate to slip instead of slip rate?
            if integrate==True:
                t=arange(0,len(stf)*stf_dt,stf_dt)
                if len(t)>len(stf):
                    t=t[0:-1]
                print(t.shape)
                print(stf.shape)
                stf_integrated=cumtrapz(stf,t,initial=0)
                stf=stf_integrated
            
            #Write stf 6 values per line
            for kstf in range(NTstf):
                if kstf==0:
                    white_space='  '
                elif (kstf+1) % 6 == 0:
                    white_space='\n'
                elif (kstf+1)==NTstf:
                    white_space='\n'
                else:
                    white_space='  '
                
                if kstf==0:
                    pre_white_space='  '
                elif (kstf) % 6 == 0:
                    pre_white_space='  '
                #elif (kstf+1)==NTstf:
                #    pre_white_space='  '
                else:
                    pre_white_space=''
                fout.write('%s%.6e%s' % (pre_white_space,stf[kstf],white_space))

    
    # And done
    print('minNTstf is: '+str(minNTstf))
    fout.close()
    




def mudpy2pecmp(rupt,outfile,thresh=0.0):
    '''
    Convert mudpy rupture to input for poroelastic PEGRN/PECMP calculation
    '''
    
    from numpy import genfromtxt,cos,deg2rad,sin
    from pyproj import Geod
    
    
    f=genfromtxt(rupt)
    #projection object
    p=Geod(ellps='WGS84')

    fout=open(outfile,'w')
    for k in range(len(f)):
        
        #read subfault quantities
        lon=f[k,1]
        lat=f[k,2]
        depth=f[k,3]
        strike=f[k,4]
        dip=f[k,5]
        ss=f[k,8]
        ds=f[k,9]
        L=f[k,10]
        W=f[k,11]
        #Get coordinates of top center
        lon_center,lat_center,az=p.fwd(lon,lat,strike-90,(W/2.)*cos(deg2rad(dip)))
        #Find the top left corner of the fault (again, godammit)
        lon_top_left,lat_top_left,az=p.fwd(lon_center,lat_center,strike-180,L/2.0)
        depth_top_left=depth-(W/2.)*sin(deg2rad(dip))/1000
        if depth_top_left<0:
            depth_top_left=0.0
        
        #Apply threshold filter
        slip=(ss**2+ds**2)**0.5
        if slip<thresh:
            ss=0
            ds=0
        
        line1='%d\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t1\t1\t0.0\n' %(k+1,lat_top_left,lon_top_left,depth_top_left,L/1000.,W/1000.,strike,dip)
        line2='\t%.2f\t%.2f\t%.2f\t%.2f\t0.0\n' %(L/2000.,W/2000.,ss,ds)
        fout.write(line1)
        fout.write(line2)
    fout.close()
    


def convolution_matrix(h):
    '''
    Return the TOeplitz-like convolution matrix for a time series
    '''
    
    from numpy import zeros,r_
    from scipy.linalg import toeplitz
    
    
    padding = zeros(h.shape[0] - 2, h.dtype)
    first_col = r_[h, padding]
    first_row = r_[h[0], padding]
    Htemp = toeplitz(first_col, first_row)
    Hfinal= Htemp[0:len(h),:]
    return Hfinal
    
    
def build_source_time_function(rise_time,dt,total_time,stf_type='triangle',zeta=0.2,dreger_falloff_rate=4,
                               scale=True,scale_value=1.0,time_offset=0,time_offset_gauss=0,quiet=False,
                               ji_fraction=0.5):
    '''
    Compute source time function for a given rise time
    '''
    from numpy import zeros,arange,where,pi,cos,sin,isnan,exp,roll
    from scipy.integrate import trapz
    
    rise_time=float(rise_time)
    #Initialize outputs
    t=arange(0,total_time+dt,dt)
    Mdot=zeros(t.shape)
    if stf_type=='triangle':
        #Triangle gradient
        m=1/(rise_time**2)
        #Upwards intercept
        b1=0
        #Downwards intercept
        b2=m*(rise_time)
        #Assign moment rate
        i=where(t<=rise_time/2)[0]
        Mdot[i]=m*t[i]+b1
        i=where(t>rise_time/2)[0]
        Mdot[i]=-m*t[i]+b2 
        i=where(t>rise_time)[0]
        Mdot[i]=0
    elif stf_type=='cosine':   #From Liu et al. (2006) BSSA, eq 7a,7b
        tau1=0.13*rise_time
        tau2=rise_time-tau1
        Cn=pi/(1.4*pi*tau1+1.2*tau1+0.3*pi*tau2)
        #Build in pieces
        i1=where(t<tau1)[0]
        i2=where((t>=tau1) & (t<2*tau1))[0]
        i3=where((t>=2*tau1) & (t<rise_time))[0]
        Mdot[i1]=Cn*(0.7-0.7*cos(t[i1]*pi/tau1)+0.6*sin(0.5*pi*t[i1]/tau1))
        Mdot[i2]=Cn*(1.0-0.7*cos(t[i2]*pi/tau1)+0.3*cos(pi*(t[i2]-tau1)/tau2))
        Mdot[i3]=Cn*(0.3+0.3*cos(pi*(t[i3]-tau1)/tau2))
    elif stf_type=='dreger':
        tau=rise_time/dreger_falloff_rate
        Mdot=(t**zeta)*exp(-t/tau)
    elif stf_type == 'gauss_prem_i_2s':  #The decay parameter for this is fixed by instaseis/syngine
    
        decay=3.5 #Hard coded in syngine, good number for rise time to actually correspond to the function width
        center_time = 3.5858 #Hard coded in syngine
        min_rise_time = 2.1 #Rise times shorter than this are not allowed
        if rise_time < min_rise_time:
            if quiet==False:
               print('... ... ...  WARNING: rise time requested is below minimum allowed of %.1fs, defaulting to minimum' % min_rise_time)
            rise_time = min_rise_time
    
        if time_offset_gauss == 0 : #regualr gaussian with no offset (for GNSS, strong motion, etc)
            time_offset_gauss = rise_time/2
            Mdot = decay/(rise_time*pi**0.5)*exp(-((decay * (t - time_offset_gauss))/rise_time)**2)
            Mdot -= Mdot[0]
            i = where(Mdot < 0)[0]
            Mdot[i] = 0
    
        else: #regular gaussian with offset (for teleseismic, etc)
            Mdot = decay/(rise_time*pi**0.5)*exp(-((decay * (t - time_offset_gauss))/rise_time)**2)
            #Here comes the tricky bit. Need to roll the slip rate function FORWARD
            # By an ammount equal to half the new rise time minus half original rise time (2.1/2)
            t_roll = rise_time/2 - min_rise_time/2
            t_roll_samples = int(t_roll/dt)
            Mdot = roll(Mdot,t_roll_samples)
    elif stf_type == 'ji':
        
        tup = rise_time*ji_fraction
        tdown = rise_time*(1-ji_fraction)
    
        #Up going cosine
        s1 = (1./(tup+tdown))*(1-cos((pi*t)/tup))
        i = where(t>tup)[0]
        s1[i] = 0
        
        #Down going cosine
        s2 = (1./(tup+tdown))*(1+cos((pi*(t-tup))/tdown))
        i = where(t<=tup)[0]
        s2[i] = 0 
        i = where(t>tup+tdown)[0]
        s2[i] = 0
        
        #add the two 
        Mdot = s1+s2    
    
    else:
        print('ERROR: unrecognized STF type '+stf_type)
        return
    #Area of STF must be equal to dt
    if scale==True:

        area=trapz(Mdot,t)
        Mdot=Mdot*(scale_value/area)
    #Check for errors
    if isnan(Mdot[0])==True:
        print('ERROR: woops, STF has nan values!')
        return
        
    #offset origin time
    t=t+time_offset
    return t,Mdot
    
    
    
def make_gnss_noise(percentile=50,dt=1,duration=600):
    
    """
    Generate synthetic GNSS noise time series data.
    
    Parameters
    ----------
    percentile : float, optional
        The percentile level to use when generating the power spectra for the noise data.
        The default value is 50.
    dt : float, optional
        The time step (in seconds) to use when generating the noise data.
        The default value is 1.
    duration : float, optional
        The total duration (in seconds) of the noise time series to generate.
        The default value is 600.
    
    Returns
    -------
    tuple
        A tuple of three arrays containing the noise time series data for the
        eastward, northward, and upward components of motion, respectively.
    
    """
    
    from mudpy.hfsims import windowed_gaussian,apply_spectrum
    
    # build white noise time series
    std=1.0 # this is a dummy parameter give it any value
    E_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    N_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    Z_noise=windowed_gaussian(duration,dt,std=std,window_type=None)
    
    # Now get the power spectra for each component of motion of noise at the 
    # specified percentile level
    f,Epsd,Npsd,Zpsd=gnss_psd(level=percentile,return_as_frequencies=True,return_as_db=False)
    
    #Covnert PSDs to amplitude spectrum
    Epsd = Epsd**0.5
    Npsd = Npsd**0.5
    Zpsd = Zpsd**0.5
    
    #apply the spectrum to the noise time series
    E_noise=apply_spectrum(E_noise,Epsd,f,dt,is_gnss=True)
    N_noise=apply_spectrum(N_noise,Npsd,f,dt,is_gnss=True)
    Z_noise=apply_spectrum(Z_noise,Zpsd,f,dt,is_gnss=True)
    
    return E_noise,N_noise,Z_noise
    
    
    
    
def read_fakequakes_hypo_time(home,project_name,rupture_name,get_Mw=False):
    '''
    Read a fakequkaes log file and extrat hypocentral time
    '''
    
    from obspy.core import UTCDateTime
    from numpy import array
    
    #old
    #rupture = rupture_name.split('.')[0]+'.'+rupture_name.split('.')[1]
    #new (Tara's fix)
    rupture = rupture_name.rsplit('.', 1)[0]
    
    log_file=home+project_name+'/output/ruptures/'+rupture+'.log'
    f=open(log_file,'r')
    while True:
        line=f.readline()
        if 'Actual magnitude' in line:
            mag=line.split('Mw')[-1]
            Mw=float(mag)
        if 'Hypocenter (lon,lat,z[km])' in line:
            s=line.split(':')[-1]
            s=s.replace('(','')
            s=s.replace(')','')
            epicenter=array(s.split(',')).astype('float')
        if 'Hypocenter time' in line:
            time_epi=line.split(' ')[-1]
            time_epi=UTCDateTime(time_epi)
            break
        if line=='':
            print('ERROR: No hypocentral time in log file')
            time_epi=''
            break
    if get_Mw==True:
        return epicenter,time_epi,Mw
    else:
        return epicenter,time_epi
    


def convert_to_grid(home,project_name,GF_list,out_file):
    '''
    Convert a bunch of single file neu's to one big grid file
    '''
    
    from numpy import genfromtxt
    
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    fout=home+'/'+project_name+'/output/forward_models/'+out_file
    f=open(fout,'w')
    
    #read offsets
    f.write('# sta,lon,lat,N(m),E(m),U(m)\n')
    for k in range(len(sta)):
        neu=genfromtxt(home+project_name+'/output/forward_models/'+sta[k]+'.static.neu')
        line='%s\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\n' % (sta[k],lonlat[k,0],lonlat[k,1],neu[0],neu[1],neu[2])
        f.write(line)
    f.close()
    
    
    
def get_Mw(rupture_file):
    '''
    Compute total moment from an inversion
    '''
    from numpy import log10,genfromtxt

   
    #Open model file
    rupt=genfromtxt(rupture_file)
    #Get slip quantities
    ss=rupt[:,8]
    ds=rupt[:,9]
    #Get total slip
    slip=(ss**2+ds**2)**0.5
    #Get rigidities and areas
    mu=rupt[:,13]
    A=rupt[:,11]*rupt[:,12]

    #Compute moments
    M0=mu*A*slip
    #Total up and copute magnitude
    M0=M0.sum()
    Mw=(2./3)*(log10(M0)-9.1)
    return M0,Mw

def gnss_noise_model():
    
    from numpy import genfromtxt
    from os import environ
    
    mudpy_aux_folder=environ['MUD']+'/src/aux/'
    
    g=genfromtxt(mudpy_aux_folder+'gnss_noise_model.txt')
    periods=g[:,0]
    psd=g[:,1:]
    
    return periods,psd

def gnss_psd(level=50,return_as_frequencies=False,return_as_db=True):
    '''
    Return PSD for all three components of motion
    '''
    
    from numpy import r_
    
    periods,psd=gnss_noise_model()
    
    if level==1:
        Epsd=psd[:,0]
        Npsd=psd[:,1]
        Zpsd=psd[:,2]
    elif level==10:
        Epsd=psd[:,3]
        Npsd=psd[:,4]
        Zpsd=psd[:,5]
    elif level==20:
        Epsd=psd[:,6]
        Npsd=psd[:,7]
        Zpsd=psd[:,8]
    elif level==30:
        Epsd=psd[:,9]
        Npsd=psd[:,10]
        Zpsd=psd[:,11]
    elif level==40:
        Epsd=psd[:,12]
        Npsd=psd[:,13]
        Zpsd=psd[:,14]
    elif level==50:
        Epsd=psd[:,15]
        Npsd=psd[:,16]
        Zpsd=psd[:,17]
    elif level==60:
        Epsd=psd[:,18]
        Npsd=psd[:,19]
        Zpsd=psd[:,20]
    elif level==70:
        Epsd=psd[:,21]
        Npsd=psd[:,22]
        Zpsd=psd[:,23]
    elif level==80:
        Epsd=psd[:,24]
        Npsd=psd[:,25]
        Zpsd=psd[:,26]
    elif level==90:
        Epsd=psd[:,27]
        Npsd=psd[:,28]
        Zpsd=psd[:,29]
        
    if return_as_frequencies==True:
        
        periods=1/periods[::-1]
        
        #reverse psds
        Epsd=Epsd[::-1]
        Npsd=Npsd[::-1]
        Zpsd=Zpsd[::-1]
        
        #add zero frequency
        periods=r_[0,periods]
        
        #add zero frequency values
        Epsd=r_[Epsd[0],Epsd]
        Npsd=r_[Npsd[0],Npsd]
        Zpsd=r_[Zpsd[0],Zpsd]
        
    if return_as_db ==False: 
        
        #Un-decibelize (not a word I know)
        Epsd=10**(Epsd/10)
        Npsd=10**(Npsd/10)
        Zpsd=10**(Zpsd/10)
        
    return periods, Epsd, Npsd, Zpsd
    
#    Ep1=





def sta_close_to_rupt(home,project_name,rupt_file,GF_list,dist_threshold,new_GF_list):
    '''
    Filt the station too far away
    IN:
        rupt_file: name of the rupture file (e.g. subduction.000036.rupt)
        GF_list: the original (large) GF_list file
        dist_threshold: keeps only the stations within this distance (unit:degree)
    OUTfile:
        new_GF_list
    '''
    import numpy as np
    import obspy
    from obspy.geodetics import kilometers2degrees
    rupt_file_path=home+project_name+'/output/ruptures/'+rupt_file
    GF_list_path=home+project_name+'/data/station_info/'+GF_list
    #load rpture file
    source=np.genfromtxt(rupt_file_path)
    rise_times=source[:,7]
    rupture_onset=source[:,12]
    #How many subfaults are non-zero?
    idxs_non_zero=np.where(rise_times>0)[0]

    #Load gflist
    STAinfo=np.genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[0,1,2],skip_header=1,dtype='S10') #higher resolution for float
    STAname=STAinfo[:,0]
    STAname=[n_staname.decode() for n_staname in STAname]
    STAlon=STAinfo[:,1]
    STAlat=STAinfo[:,2]
    STAlon=np.array([float(i) for i in STAlon])
    STAlat=np.array([float(i) for i in STAlat])
    OUT1=open(home+project_name+'/data/station_info/'+new_GF_list,'w')
    OUT1.write('#\n') #header
    #loop through all stations
    for nsta in range(len(STAname)):
        for idx_non_zero in idxs_non_zero:
            sub_lon=source[idx_non_zero,1]
            sub_lat=source[idx_non_zero,2]
            dist=obspy.geodetics.locations2degrees(lat1=sub_lat,long1=sub_lon,lat2=STAlat[nsta],long2=STAlon[nsta])
            #dist,az,baz=obspy.geodetics.base.gps2dist_azimuth(lat1=sub_lat,lon1=sub_lon,lat2=STAlat[nsta],lon2=STAlon[nsta])
            if dist<dist_threshold:
                OUT1.write('%s %10.6f %10.6f 0 1 0 0 0 0 /foo/bar /foo/bar /foo/bar /foo/bar /foo/bar 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n'%(STAname[nsta],STAlon[nsta],STAlat[nsta]))
                break
    OUT1.close()

