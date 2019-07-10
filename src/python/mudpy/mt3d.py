'''
Routines for 3d or volumetric moment tensor inversion

D. Melgar
Apr 4th 2017
'''



#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,source_name,model_name,
        green_flag,synth_flag,dk,pmin,pmax,kmax,hot_start,ncpus,forceMT,mt):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,where
    from os import remove
    from gc import collect
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='U')
    GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
    
    #Parameters not used but needed in function call
    dt=1.0
    NFFT=64
    
    # GFs can be computed all at the same time
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    if green_flag==1:
        #decide what GF computation is required for this station
        i=where(GF[:,2]==1)[0]
        
        if len(i)>0: #Static offset
            print('Static GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_green(home,project_name,station_file,source_name,model_name,dt,NFFT,
                            hot_start,dk,pmin,pmax,kmax,ncpus,insar=False)
        
        i=where(GF[:,6]==1)[0]
        if len(i)>0: #InSAR LOS
            print('InSAR GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_green(home,project_name,station_file,source_name,model_name,dt,NFFT,
                            hot_start,dk,pmin,pmax,kmax,ncpus,insar=True)

            collect()   
    
    # Synthetics
    if synth_flag==1:

        station_file='temp.sta'
        #Decide which synthetics are required
        
        i=where(GF[:,2]==1)[0]
        if len(i)>0: #Static offset
            print('Static synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus,forceMT,mt,insar=False)

        #Decide which synthetics are required
        i=where(GF[:,6]==1)[0]
        if len(i)>0: # InSAR LOS
            print('InSAR synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus,forceMT,mt,insar=True)





def make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,
            hot_start,dk,pmin,pmax,kmax,ncpus,insar=False):
    '''
    This routine set's up the computation of GFs for each subfault to all stations.
    The GFs are impulse sources, they don't yet depend on strike and dip.
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        dt: Desired sampling itnerval for waveforms
        NFFT: No. of samples requested in waveform (must be power of 2)
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all

        
    OUT:
        Nothing
    '''
    from numpy import loadtxt,arange,savetxt
    from os import path,makedirs,environ
    from shlex import split
    import subprocess
    
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    #Load source model for station-event distance computations
    source=loadtxt(fault_file,ndmin=2)

    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    #Make mpi system call
    print("MPI: Starting GFs computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'

    #Parameters that are not used but neede by function call
    static=1
    tsunami=False
    
    #Make mpi call
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_green '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(dt)+' '+str(NFFT)+' '+str(static)+' '+str(dk)+' '+str(pmin)+' '+str(pmax)+' '+str(kmax)+' '+str(tsunami)+' '+str(insar)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
    
    
    
    
#Now make synthetics for source/station pairs
def make_parallel_synthetics(home,project_name,station_file,source_name,model_name,
                    hot_start,ncpus,forceMT,mt,insar=False):
    '''
    This routine will take the impulse response (GFs) and pass it into the routine that will
    convovle them with the source time function according to each subfaults strike and dip.
    The result fo this computation is a time series dubbed a "synthetic"
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to be displacement
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        Nothing
    '''
    from numpy import arange,savetxt
    import datetime
    from numpy import loadtxt
    import subprocess
    from shlex import split
    from os import environ
    
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+source_name
    synth_path=home+project_name+'/GFs/static/'
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    
    if forceMT==False:
        mt=[0,0,0,0,0,0]
    
    #Make mpi system call
    print("MPI: Starting synthetics computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_synthetics_mt3d '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(forceMT)+' '+str(mt[0])+' '+str(mt[1])+' '+str(mt[2])+' '+str(mt[3])+' '+str(mt[4])+' '+str(mt[5])+' '+str(insar)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
    
    #Now merge files (loop over MT components
    #Make list of files to be merged
    filenames_Mxx=[]
    filenames_Mxy=[]
    filenames_Mxz=[]
    filenames_Myy=[]
    filenames_Myz=[]
    filenames_Mzz=[]
    
    #List of file names to cocnatenate
    for kfile in range(ncpus):
        filenames_Mxx.append(synth_path+'tmp_Mxx_process'+str(int(kfile)))
        filenames_Mxy.append(synth_path+'tmp_Mxy_process'+str(int(kfile)))
        filenames_Mxz.append(synth_path+'tmp_Mxz_process'+str(int(kfile)))
        filenames_Myy.append(synth_path+'tmp_Myy_process'+str(int(kfile)))
        filenames_Myz.append(synth_path+'tmp_Myz_process'+str(int(kfile)))
        filenames_Mzz.append(synth_path+'tmp_Mzz_process'+str(int(kfile)))

    
    #Make final output files for GPS
    if forceMT==True:
        if insar==False:
            with open(synth_path+'_Mxx.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Mxx:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
    
        else: #Output for Insar
            with open(synth_path+'_Mxx.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Mxx:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            

        
    else: #Stuff to do when you want all 6 MT components
        if insar==False:
            with open(synth_path+'_Mxx.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Mxx:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mxy.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n') 
                for fname in filenames_Mxy:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mxz.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Mxz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Myy.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Myy:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Myz.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Myz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mzz.neu', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),n(m),e(m),u(m)\n')
                for fname in filenames_Mzz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
        
    
        else: #Output for Insar
            with open(synth_path+'_Mxx.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Mxx:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mxy.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n') 
                for fname in filenames_Mxy:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mxz.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Mxz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Myy.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Myy:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Myz.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Myz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                            
            with open(synth_path+'_Mzz.los', 'w') as outfile:
                outfile.write('# sta_name,src_num,lon_sta,lat_sta,lon_src,lat_src,z_src(km),los(m)\n')
                for fname in filenames_Mzz:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
                        
    

def getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name,forceMT):
    '''
    Build the G matrix
    '''
    from numpy import arange,genfromtxt,where,loadtxt,array,c_,concatenate,save,load,size,tile,expand_dims
    from os import remove
    from os.path import split
    from mudpy import inverse as inv
    
    G_name=home+project_name+'/GFs/matrices/'+G_name
    if G_from_file==True: #load from file
        print('Loading G from file '+G_name)
        G=load(G_name+'.npy')
    else: #assemble G one data type at a time
        print('Assembling G from synthetic computations...')
        #Read in GFlist and decide what to compute
        gf_file=home+project_name+'/data/station_info/'+GF_list
        mini_station=home+project_name+'/data/station_info/tempG.sta'
        stations=genfromtxt(gf_file,usecols=0,dtype='U')
        GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],dtype='f8')
        GFfiles=genfromtxt(gf_file,usecols=[8,9,10,11,12],dtype='S')
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
                inv.mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,0])
                gftype='static'
                Gstatic= makeG(home,project_name,source_name,model_name,split(mini_station)[1],gftype,forceMT)
                remove(mini_station) #Cleanup  

        #static field GFs
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
                inv.mini_station_file(mini_station,stations[i],GF[i,0],GF[i,1],GFfiles[i,0])
                gftype='insar'
                Ginsar= makeG(home,project_name,source_name,model_name,split(mini_station)[1],gftype,forceMT)
                remove(mini_station) #Cleanup  
        
        #Done, concat, save to file, release the memory
        print(Gstatic.shape)
        print(Ginsar.shape)
        G=concatenate([g for g in [Gstatic,Ginsar] if g.size > 0])
        print('Saving GF matrix to '+G_name+' this might take just a second...')
        save(G_name,G)
    return G



def makeG(home,project_name,source_name,model_name,station_file,gftype,forceMT):
    '''
    This routine is called from getG and will assemble the GFs from available synthetics
    depending on data type requested 
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        model_name: Name of velocity structure file
        station_file: File with coordinates of stations and data types
        gftype: ='static' if assembling static field GFs
       
    OUT:
        G: Partially assembled GF with all synthetics from a particular data type
    '''
    
    from numpy import genfromtxt,loadtxt,zeros,array,inf,where,unique,argsort,arange

    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+source_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="U",usecols=0)
    datafiles=genfromtxt(station_file,dtype="S",usecols=3)
    syn_path=home+project_name+'/GFs/static/'
    
    #Deal with one station issue
    try:
        Nsta=len(staname)
    except:
        Nsta=1
        staname=array([staname])
        datafiles=array([datafiles])
    
    #Initalize G for faster assignments
    if forceMT==True:
        Ncomponents=1
    else:
        Ncomponents=6
        
    if gftype.lower()=='static': #Initialize output matrix
        G=zeros((Nsta*3,Nfaults*Ncomponents))
    elif gftype.lower()=='insar': #Initialize output matrix
        G=zeros((Nsta,Nfaults*Ncomponents))
        
    if gftype.lower()=='static': #Make matrix of static GFs
        #Load the synthetics files
        Mxx=genfromtxt(syn_path+'_Mxx.neu')
        Mxy=genfromtxt(syn_path+'_Mxy.neu')
        Mxz=genfromtxt(syn_path+'_Mxz.neu')
        Myy=genfromtxt(syn_path+'_Myy.neu')
        Myz=genfromtxt(syn_path+'_Myz.neu')
        Mzz=genfromtxt(syn_path+'_Mzz.neu')
        Mxx_sta=genfromtxt(syn_path+'_Mxx.neu',usecols=0,dtype='S')
        Mxy_sta=genfromtxt(syn_path+'_Mxy.neu',usecols=0,dtype='S')
        Mxz_sta=genfromtxt(syn_path+'_Mxz.neu',usecols=0,dtype='S')
        Myy_sta=genfromtxt(syn_path+'_Myy.neu',usecols=0,dtype='S')
        Myz_sta=genfromtxt(syn_path+'_Myz.neu',usecols=0,dtype='S')
        Mzz_sta=genfromtxt(syn_path+'_Mzz.neu',usecols=0,dtype='S')
        
        for ksta in range(Nsta):
            
            #Loop over stations
            if ksta%10==0:
                print('... working on station '+str(ksta)+' of '+str(Nsta))
            
            #Initalize output variable
            if forceMT==False:
                Gtemp=zeros([3,Nfaults*6])
            else:
                Gtemp=zeros([3,Nfaults])
                
            #what am I working on
            sta=staname[ksta]
            
            #Find the position ine ach file of allsubfaults for this station
            i=where(Mxx_sta==sta)[0]
            isort=argsort(Mxx[i,1])
            neu_xx=Mxx[i[isort],7:10]
            
            i=where(Mxy_sta==sta)[0]
            isort=argsort(Mxy[i,1])
            neu_xy=Mxy[i[isort],7:10]
            
            i=where(Mxz_sta==sta)[0]
            isort=argsort(Mxz[i,1])
            neu_xz=Mxz[i[isort],7:10]
            
            i=where(Myy_sta==sta)[0]
            isort=argsort(Myy[i,1])
            neu_yy=Myy[i[isort],7:10]
            
            i=where(Myz_sta==sta)[0]
            isort=argsort(Myz[i,1])
            neu_yz=Myz[i[isort],7:10]
            
            i=where(Mzz_sta==sta)[0]
            isort=argsort(Mzz[i,1])
            neu_zz=Mzz[i[isort],7:10]
            
            ##Needs work
            ##Place into G matrix
            #if forceMT==False:
            #    #North
            #    Gtemp[0,6*kfault]=neu_xx[0,0]   ; Gtemp[0,6*kfault+1]=neu_xy[0,0] ; Gtemp[0,6*kfault+2]=neu_xz[0,0]
            #    Gtemp[0,6*kfault+3]=neu_yy[0,0] ; Gtemp[0,6*kfault+4]=neu_yz[0,0] ; Gtemp[0,6*kfault+5]=neu_zz[0,0]
            #    
            #    #East
            #    Gtemp[1,6*kfault]=neu_xx[0,1]   ; Gtemp[1,6*kfault+1]=neu_xy[0,1] ; Gtemp[1,6*kfault+2]=neu_xz[0,1]
            #    Gtemp[1,6*kfault+3]=neu_yy[0,1] ; Gtemp[1,6*kfault+4]=neu_yz[0,1] ; Gtemp[1,6*kfault+5]=neu_zz[0,1]
            #    
            #    #Up
            #    Gtemp[2,6*kfault]=neu_xx[0,2]   ; Gtemp[2,6*kfault+1]=neu_xy[0,2] ; Gtemp[2,6*kfault+2]=neu_xz[0,2]
            #    Gtemp[2,6*kfault+3]=neu_yy[0,2] ; Gtemp[2,6*kfault+4]=neu_yz[0,2] ; Gtemp[2,6*kfault+5]=neu_zz[0,2]
            #else:
            #    #North
            #    Gtemp[0,kfault]=neu_xx[0,0]
            #    #East
            #    Gtemp[1,kfault]=neu_xx[0,1]                 
            #    #Up
            #    Gtemp[2,kfault]=neu_xx[0,2]
            #    
            ##Append to output matrix
            #G[ksta*3:ksta*3+3,:]=Gtemp   
        return G   
    
    elif gftype.lower()=='insar': #Make matrix of insar LOS GFs
        #Load the synthetics files
        print('... loading Mxx')
        Mxx=genfromtxt(syn_path+'_Mxx.los')
        Mxx_sta=genfromtxt(syn_path+'_Mxx.los',usecols=0,dtype='S')
        if forceMT==False:
            print('... loading Mxy')
            Mxy=genfromtxt(syn_path+'_Mxy.los')
            Mxy_sta=genfromtxt(syn_path+'_Mxy.los',usecols=0,dtype='S')
            print('... loading Mxz')
            Mxz=genfromtxt(syn_path+'_Mxz.los')
            Mxz_sta=genfromtxt(syn_path+'_Mxz.los',usecols=0,dtype='S')
            print('... loading Myy')
            Myy=genfromtxt(syn_path+'_Myy.los')
            Myy_sta=genfromtxt(syn_path+'_Myy.los',usecols=0,dtype='S')
            print('... loading Myz')
            Myz=genfromtxt(syn_path+'_Myz.los')
            Myz_sta=genfromtxt(syn_path+'_Myz.los',usecols=0,dtype='S')
            print('... loading Mzz')
            Mzz=genfromtxt(syn_path+'_Mzz.los')
            Mzz_sta=genfromtxt(syn_path+'_Mzz.los',usecols=0,dtype='S')     
        
        
        for ksta in range(Nsta):
            #Loop over stations
            if ksta%10==0:
                print('... working on station '+str(ksta)+' of '+str(Nsta))
            
            #Initalize output variable
            if forceMT==False:
                Gtemp=zeros([1,Nfaults*6])
            else:
                Gtemp=zeros([1,Nfaults])
                
            #what am I working on
            sta=staname[ksta]
            
            #Find the position ine ach file of allsubfaults for this station
            i=where(Mxx_sta==sta)[0]
            isort=argsort(Mxx[i,1])
            neu_xx=Mxx[i[isort],7]
            
            i=where(Mxy_sta==sta)[0]
            isort=argsort(Mxy[i,1])
            neu_xy=Mxy[i[isort],7]
            
            i=where(Mxz_sta==sta)[0]
            isort=argsort(Mxz[i,1])
            neu_xz=Mxz[i[isort],7]
            
            i=where(Myy_sta==sta)[0]
            isort=argsort(Myy[i,1])
            neu_yy=Myy[i[isort],7]
            
            i=where(Myz_sta==sta)[0]
            isort=argsort(Myz[i,1])
            neu_yz=Myz[i[isort],7]
            
            i=where(Mzz_sta==sta)[0]
            isort=argsort(Mzz[i,1])
            neu_zz=Mzz[i[isort],7]
                    
            #Place into G matrix
            # LOS
            if forceMT==False:
                ixx=arange(0,Nfaults*6,6)
                ixy=arange(1,Nfaults*6,6)
                ixz=arange(2,Nfaults*6,6)
                iyy=arange(3,Nfaults*6,6)
                iyz=arange(4,Nfaults*6,6)
                izz=arange(5,Nfaults*6,6)
                Gtemp[0,ixx]=neu_xx   ; Gtemp[0,ixy]=neu_xy ; Gtemp[0,ixz]=neu_xz
                Gtemp[0,iyy]=neu_yy ; Gtemp[0,iyz]=neu_yz ; Gtemp[0,izz]=neu_zz
            else:
                Gtemp[0,kfault]=neu_xx

                
            #Append to output matrix
            G[ksta,:]=Gtemp   
        return G   


def get_moment(home,project_name,source_name,sol,forceMT):
    '''
    Compute total moment from an inversion
    '''
    from numpy import log10,genfromtxt,loadtxt,arange,zeros,array
    from mudpy.forward import get_mu
    from numpy.linalg import eig
   
    unitMw=5.0
    unitM0=10**(unitMw*1.5+9.1) #Keep it in N-m
   
    if forceMT==True: 
        M0=sol*unitM0  
    
    else:
        #Open model file
        f=genfromtxt(home+project_name+'/data/model_info/'+source_name)
    
        #Get mt components
        ixx=6*arange(len(sol)/6)
        ixy=6*arange(len(sol)/6)+1
        ixz=6*arange(len(sol)/6)+2
        iyy=6*arange(len(sol)/6)+3
        iyz=6*arange(len(sol)/6)+4
        izz=6*arange(len(sol)/6)+5
        
        mt_xx=sol[ixx,0]
        mt_xy=sol[ixy,0]
        mt_xz=sol[ixz,0]
        mt_yy=sol[iyy,0]
        mt_yz=sol[iyz,0]
        mt_zz=sol[izz,0]
        
        #Get total moment
        M0=zeros(len(mt_xx))
        for k in range(len(mt_xx)):
            M=array([[mt_xx[k],mt_xy[k],mt_xz[k]],[mt_xy[k],mt_yy[k],mt_yz[k]],[mt_xz[k],mt_yz[k],mt_zz[k]]])
            M=M*unitM0
            eigVal,eigVec=eig(M)
            M0[k]=(0.5*sum(eigVal**2))**0.5
        
    
    #Total up and copute magnitude
    M0total=M0.sum()
    Mw=(2./3)*(log10(M0total)-9.1)
    
    return M0total,M0,Mw


def write_log(home,project_name,run_name,k,lambda_spatial,
        L2,Lm,VR,Mo,Mw,velmod,source,g_name,gflist):
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
    f.write('Source model: '+source+'\n')
    f.write('G name: '+g_name+'\n')
    f.write('GF list: '+gflist+'\n')
    f.write('lambda_spatial = '+repr(lambda_spatial)+'\n')
    f.write('L2 = '+repr(L2)+'\n')
    f.write('VR static(%) = '+repr(VR[0])+'\n')
    f.write('VR InSAR LOS(%) = '+repr(VR[4])+'\n')
    f.write('Lm = '+repr(Lm)+'\n')
    f.write('M0(N-m) = '+repr(Mo)+'\n')
    f.write('Mw = '+repr(Mw)+'\n')
    f.close()



def write_model(home,project_name,run_name,source_name,model_name,sol,num,forceMT,mt):
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
    
    from numpy import genfromtxt,loadtxt,arange,zeros,c_,savetxt,r_,array
    from mudpy.forward import get_mu
    from numpy.linalg import eig
    from mudpy import analysis
   
    #Open model file
    f=genfromtxt(home+project_name+'/data/model_info/'+source_name)
    
    unitMw=5.0
    unitM0=10**(unitMw*1.5+9.1) #Keep it in N-m
    
    #Get MT components and MT related quantities (moment, strike dip, etc)-
    ixx=6*arange(len(sol)/6)
    ixy=6*arange(len(sol)/6)+1
    ixz=6*arange(len(sol)/6)+2
    iyy=6*arange(len(sol)/6)+3
    iyz=6*arange(len(sol)/6)+4
    izz=6*arange(len(sol)/6)+5
    
    if forceMT==True:
        mt_xx=mt[0]*sol*unitM0
        mt_xy=mt[1]*sol*unitM0
        mt_xz=mt[2]*sol*unitM0
        mt_yy=mt[3]*sol*unitM0
        mt_yz=mt[4]*sol*unitM0
        mt_zz=mt[5]*sol*unitM0
    
    else: 
        mt_xx=sol[ixx,0]*unitM0
        mt_xy=sol[ixy,0]*unitM0
        mt_xz=sol[ixz,0]*unitM0
        mt_yy=sol[iyy,0]*unitM0
        mt_yz=sol[iyz,0]*unitM0
        mt_zz=sol[izz,0]*unitM0
    
    #Get total moment
    M0=zeros(len(mt_xx))
    strike1=zeros(len(mt_xx))
    dip1=zeros(len(mt_xx))
    rake1=zeros(len(mt_xx))
    strike2=zeros(len(mt_xx))
    dip2=zeros(len(mt_xx))
    rake2=zeros(len(mt_xx))
    for k in range(len(mt_xx)):
        mt=analysis.MT(mt_xx[k],mt_xy[k],mt_xz[k],mt_yy[k],mt_yz[k],mt_zz[k],1,1,1,mt_style='xyz')
        mt.get_nodal_planes()
        NP1=mt.nodal_plane1
        NP2=mt.nodal_plane2
        M0[k]=mt.moment()
        strike1[k]=NP1[0]
        dip1[k]=NP1[1]
        rake1[k]=NP1[2]
        strike2[k]=NP2[0]
        dip2[k]=NP2[1]
        rake2[k]=NP2[2]
    

    #Prepare for output
    out=c_[f,mt_xx,mt_xy,mt_xz,mt_yy,mt_yz,mt_zz,strike1,dip1,rake1,strike2,dip2,rake2,M0]
    outdir=home+project_name+'/output/inverse_models/models/'+run_name+'.'+str(num).rjust(4,'0')+'.inv'
    #CHANGE this to rupture definition as #No  x            y        z(km)      str     dip      rake       rise    dura     slip    ss_len  ds_len rupt_time
    fmtout='%6i\t%11.4f\t%11.4f\t%11.4f\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%13.4e'
    print('... writing model results to file '+outdir)
    savetxt(outdir,out,fmtout,header='No,lon,lat,z(km),Mxx(Nm),Mxy(Nm),Mxz(Nm),Myy(Nm),Myz(Nm),Mzz(Nm),strike1,dip1,rake1,strike2,dip2,rake2,M0(Nm)')



def run_inversion(home,project_name,run_name,source_name,model_name,GF_list,G_from_file,G_name,
        reg_spatial,nsources,solver,forceMT,mt,weight=False,Ltype=0):
    '''
    Assemble G and d, determine smoothing and run the inversion
    '''
    from mudpy import inverse as inv
    from mudpy.forward import get_mu_and_area
    from numpy import zeros,dot,array,squeeze,expand_dims,empty,tile,eye,ones,arange,load
    from numpy.linalg import lstsq
    from scipy.sparse import csr_matrix as sparse
    from scipy.optimize import nnls
    from datetime import datetime
    import gc
    
    

    t1=datetime.now()
    #Get data vector
    d=inv.getdata(home,project_name,GF_list,decimate=None,bandpass=None)

    #Get GFs
    G=getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name,forceMT)
    
    #Get data weights
    if weight==True:
        print('Applying data weights')
        w=inv.get_data_weights(home,project_name,GF_list,d,None)
        W=empty(G.shape)
        W=tile(w,(G.shape[1],1)).T
        WG=empty(G.shape)
        WG=W*G
        wd=w*d.squeeze()
        wd=expand_dims(wd,axis=1)
        #Clear up extraneous variables
        W=None
        w=None
        #Define inversion quantities
        x=WG.transpose().dot(wd)
        print('Computing G\'G')
        K=(WG.T).dot(WG)
        #And cleanup
        W=None
        WG=None
        wd=None
    else:
        #Define inversion quantities if no weightd
        x=G.transpose().dot(d)
        print('Computing G\'G')
        K=(G.T).dot(G)
    
    #Cleanup
    #G=None 
    
    #Get regularization matrices (set to 0 matrix if not needed)
    if forceMT==False:
        Ls=eye(nsources*6) 
    else:
        Ls=eye(nsources)
    print('Nsources: '+str(nsources))
    Ninversion=len(reg_spatial)

    #Make L's sparse
    Ls=sparse(Ls)
    #Get regularization tranposes for ABIC
    LsLs=Ls.transpose().dot(Ls)
    #inflate
    Ls=Ls.todense()
    LsLs=LsLs.todense()

    #off we go
    dt=datetime.now()-t1
    print('Preprocessing wall time was '+str(dt))
    print('\n--- RUNNING INVERSIONS ---\n')
    ttotal=datetime.now()
    kout=0
    
    for ks in range(len(reg_spatial)):
        
        t1=datetime.now()
        lambda_spatial=reg_spatial[ks]
        print('Running inversion '+str(kout+1)+' of '+str(Ninversion)+' at regularization levels: ls ='+repr(lambda_spatial))
        
        Kinv=K+(lambda_spatial**2)*LsLs
        
        if solver=='lstsq':
            sol,res,rank,s=lstsq(Kinv,x)
        elif solver=='nnls':
            sol,res=nnls(Kinv,squeeze(x.T))

        #Compute synthetics
        ds=dot(G,sol)
        
        #Get stats
        L2,Lmodel=inv.get_stats(Kinv,sol,x)
        VR=inv.get_VR(home,project_name,GF_list,sol,d,ds,None)

        #Get moment
        M0total,M0,Mw=get_moment(home,project_name,source_name,sol,forceMT)

        #Write log
        write_log(home,project_name,run_name,kout,lambda_spatial,
            L2,Lmodel,VR,M0total,Mw,model_name,source_name,G_name,GF_list)
        
        #Write output to file
        inv.write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,kout,None)
        write_model(home,project_name,run_name,source_name,model_name,sol,kout,forceMT,mt)
        
        kout+=1
        dt1=datetime.now()-t1
        dt2=datetime.now()-ttotal
        print('... inversion wall time was '+str(dt1)+', total wall time elapsed is '+str(dt2))