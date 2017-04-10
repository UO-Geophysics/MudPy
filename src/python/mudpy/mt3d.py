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
                            hot_start,dk,pmin,pmax,kmax,ncpus)
        
        i=where(GF[:,6]==1)[0]
        if len(i)>0: #InSAR LOS
            print 'InSAR GFs requested...'
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            make_parallel_green(home,project_name,station_file,source_name,model_name,dt,NFFT,
                            hot_start,dk,pmin,pmax,kmax,ncpus)

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
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus)

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
            make_parallel_synthetics(home,project_name,station_file,source_name,model_name,hot_start,ncpus)





def make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,
            hot_start,dk,pmin,pmax,kmax,ncpus):
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
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_green '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(dt)+' '+str(NFFT)+' '+str(static)+' '+str(dk)+' '+str(pmin)+' '+str(pmax)+' '+str(kmax)+' '+str(tsunami)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
    
    
    
    
#Now make synthetics for source/station pairs
def make_parallel_synthetics(home,project_name,station_file,source_name,model_name,
                    hot_start,ncpus):
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
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_synthetics_mt3d '+home+' '+project_name+' '+station_file+' '+model_name
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()