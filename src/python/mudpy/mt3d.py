'''
Routines for 3d or volumetric moment tensor inversion

D. Melgar
Apr 4th 2017
'''



#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,source_name,model_name,
        green_flag,synth_flag,dk,pmin,pmax,kmax,hot_start,ncpus):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,where
    from os import remove
    from gc import collect
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='S6')
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
            print 'Static GFs requested...'
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_green(home,project_name,station_file,source_name,model_name,dt,NFFT,
                            hot_start,dk,pmin,pmax,kmax,ncpus,insar=False)
        
        i=where(GF[:,6]==1)[0]
        if len(i)>0: #InSAR LOS
            print 'InSAR GFs requested...'
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
            print 'Static synthetics requested'
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus,insar=False)

        #Decide which synthetics are required
        i=where(GF[:,6]==1)[0]
        if len(i)>0: # InSAR LOS
            print 'InSAR synthetics requested'
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus,insar=True)





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
    from string import rjust
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
    print "MPI: Starting GFs computation on", ncpus, "CPUs\n"
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
                    hot_start,ncpus,insar=False):
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
    
    #Make mpi system call
    print "MPI: Starting synthetics computation on", ncpus, "CPUs\n"
    mud_source=environ['MUD']+'/src/python/mudpy/'
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_synthetics_mt3d '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(insar)
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
                        
    

def getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name):
    '''
    Build the G matrix
    '''
    from numpy import arange,genfromtxt,where,loadtxt,array,c_,concatenate,save,load,size,tile,expand_dims
    from os import remove
    from os.path import split
    from mudpy import inverse as inv
    
    G_name=home+project_name+'/GFs/matrices/'+G_name
    if G_from_file==True: #load from file
        print 'Loading G from file '+G_name
        G=load(G_name)
    else: #assemble G one data type at a time
        print 'Assembling G from synthetic computations...'
        #Read in GFlist and decide what to compute
        gf_file=home+project_name+'/data/station_info/'+GF_list
        mini_station=home+project_name+'/data/station_info/tempG.sta'
        stations=genfromtxt(gf_file,usecols=0,dtype='S6')
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
                tdelay=0
                Gstatic= makeG(home,project_name,source_name,model_name,split(mini_station)[1],gftype)
                remove(mini_station) #Cleanup  





def makeG(home,project_name,source_name,model_name,station_file,gftype):
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
    
    from numpy import genfromtxt,loadtxt,zeros,array,inf
    from string import rjust

    
    #Load fault model
    source=loadtxt(home+project_name+'/data/model_info/'+source_name,ndmin=2)
    Nfaults=source.shape[0] #Number of subfaults
    #Load station info
    station_file=home+project_name+'/data/station_info/'+station_file
    staname=genfromtxt(station_file,dtype="S6",usecols=0)
    datafiles=genfromtxt(station_file,dtype="S",usecols=3)
    
    #Deal with one station issue
    try:
        Nsta=len(staname)
    except:
        Nsta=1
        staname=array([staname])
        datafiles=array([datafiles])
    insert_position=0
    
    #Initalize G for faster assignments
    if gftype.lower()=='static': #Initialize output matrix
        G=zeros((Nsta*3,Nfaults*6))
    elif gftype.lower()=='insar': #Initialize output matrix
        G=zeros((Nsta,Nfaults*6))
    else:
        pass #For disp or vel waveforms G is initalized below
        
    if gftype.lower()=='static': #Make matrix of static GFs
        for ksta in range(Nsta):
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
        return G   
    if gftype.lower()=='insar': #Make matrix of insar LOS GFs
        for ksta in range(Nsta):
            print 'Assembling static GFs for station '+staname[ksta]
            #Initalize output variable
            Gtemp=zeros([1,Nfaults*2])
            #Where's the data
            syn_path=home+project_name+'/GFs/static/'
            #Data path, need this to find LOS vector
            los_path=home+project_name+'/data/statics/'
            #Read los vector for this subfault
            los=genfromtxt(los_path+staname[ksta]+'.los')
            los=los[1:]
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
                # Dot product of GFs and los vector
                los_ss=los.dot(array([nss,ess,zss]))
                los_ds=los.dot(array([nds,eds,zds]))
                #Place into G matrix
                Gtemp[0,2*kfault]=los_ss   ; Gtemp[0,2*kfault+1]=los_ds  
                #Append to G
            #Append to output matrix
            G[ksta*1,:]=Gtemp   
        return G 






def run_inversion(home,project_name,run_name,source_name,model_name,GF_list,G_from_file,G_name,
        reg_spatial,reg_temporal,nfaults,solver,weight=False,Ltype=2):
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
    G=getG(home,project_name,source_name,model_name,GF_list,G_from_file,G_name)