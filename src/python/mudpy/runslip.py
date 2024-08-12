'''
Diego Melgar, 01/2014
Runtime file for forward modeling and inverse kinematic slip inversions
'''


#Initalize project folders
def init(home,project_name):
    '''
    Initalizes file structure for a new problem
    
    IN:
        home: What dir will you be working from
        project_name: What name will you give this problem
        
    OUT:
        Nothing
    '''
    from shutil import rmtree,copy
    from os import path,makedirs,environ
    clob='y'
    proj_dir=home+project_name+'/'
    if path.exists(proj_dir):  #Path exists, clobber?
        clob=input('Project directory exists, clobber (y/n)?')
        if clob =='y' or clob == 'Y': #Clobber baby
            clob=input('This will delete everything in this project directory, so, take a minute, think about it: clobber (y/n)?')
            if clob == 'y' or clob == 'Y':
                rmtree(proj_dir)
            else: #Leave direcory alone
                print('Phew, almost shot yourself in the foot there didn\'t you?')
        else: #Leave directory alone
            print('Phew, almost shot yourself in the foot there didn\'t you?')
    if clob == 'y' or clob == 'Y':
        makedirs(proj_dir)
        #And make the subdirectories
        makedirs(proj_dir+'GFs')
        makedirs(proj_dir+'GFs/static')
        makedirs(proj_dir+'GFs/dynamic')
        makedirs(proj_dir+'GFs/matrices')
        makedirs(proj_dir+'GFs/tsunami')
        makedirs(proj_dir+'GFs/STFs')
        makedirs(proj_dir+'data/waveforms')
        makedirs(proj_dir+'data/statics')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'scripts')
        makedirs(proj_dir+'forward_models')
        makedirs(proj_dir+'output/inverse_models')
        makedirs(proj_dir+'output/inverse_models/statics')
        makedirs(proj_dir+'output/inverse_models/waveforms')
        makedirs(proj_dir+'output/inverse_models/models')
        makedirs(proj_dir+'output/forward_models')
        makedirs(proj_dir+'logs')
        makedirs(proj_dir+'analysis')
        makedirs(proj_dir+'analysis/frequency')
        #Copy templates into appropriate files
        mudpy=environ['MUD']+'/run/'
        copy(mudpy+'template.fault',proj_dir+'data/model_info/')
        copy(mudpy+'template.gflist',proj_dir+'data/station_info/')
        copy(mudpy+'template.sta',proj_dir+'data/station_info/')
        copy(mudpy+'template.mod',proj_dir+'structure/')


#Extract fault geometry from rupture file
def rupt2fault(home,project_name,rupture_name):
    '''
    Make fault file from user provided forward model rupture file
    '''
    from numpy import loadtxt,savetxt,c_
    
    print('Assembling fault file from rupture file')
    rupt=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    fault=c_[rupt[:,0],rupt[:,1],rupt[:,2],rupt[:,3],rupt[:,4],rupt[:,5],rupt[:,6],rupt[:,7],rupt[:,10],rupt[:,11]]
    savetxt(home+project_name+'/data/model_info/'+rupture_name.split('.')[0]+'.fault', \
            fault,fmt='%i\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f')





# Run green functions          
def make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
            hot_start,dk,pmin,pmax,kmax,okada=False):
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
    import time,glob
    from mudpy import green
    from numpy import loadtxt
    from shutil import rmtree,copy
    from os import chdir,path,makedirs,remove
    import datetime
    import gc
    
    if okada==False:
        tic=time.time()
        model_path=home+project_name+'/structure/'
        green_path=home+project_name+'/GFs/'
        station_file=home+project_name+'/data/station_info/'+station_file 
        fault_file=home+project_name+'/data/model_info/'+fault_name  
        logpath=home+project_name+'/logs/'
        #log time
        now=datetime.datetime.now()
        now=now.strftime('%b-%d-%H%M')
        chdir(model_path)
        #Load source model for station-event distance computations
        source=loadtxt(fault_file,ndmin=2)
        for k in range(hot_start,source.shape[0]):
            #Run comptuation for 1 subfault
            log=green.run_green(source[k,:],station_file,model_name,dt,NFFT,static,dk,pmin,pmax,kmax)
            #Write log
            f=open(logpath+'make_green.'+now+'.log','a')    
            f.write(log)
            f.close()
            #Move to correct directory
            strdepth='%.4f' % source[k,3]
            subfault=str(int(source[k,0])).rjust(4,'0') #subfault=rjust(str(int(source[k,0])),4,'0')
            if static==0 and tsunami==False:
                #Move results to dynamic GF dir
                dirs=glob.glob('*.mod_'+strdepth)
                #Where am I writting this junk too?
                outgreen=green_path+'/dynamic/'+path.split(dirs[0])[1]+'.sub'+subfault
                #Check if GF subdir already exists
                if path.exists(outgreen)==False:
                    #It doesn't, make it, don't be lazy
                    makedirs(outgreen)
                #Now copy GFs in, this will OVERWRITE EXISTING GFs of the same name
                flist=glob.glob(dirs[0]+'/*')
                for k in range(len(flist)):
                    copy(flist[k],outgreen)
                #Cleanup
                rmtree(dirs[0])
                gc.collect()
            elif static==0 and tsunami==True: #Tsunami GFs
                #Move results to tsunami GF dir
                dirs=glob.glob('*.mod_'+strdepth)
                #Where am I writting this junk too?
                outgreen=green_path+'/tsunami/'+path.split(dirs[0])[1]+'.sub'+subfault
                #Check if GF subdir already exists
                if path.exists(outgreen)==False:
                    #It doesn't, make it, don't be lazy
                    makedirs(outgreen)
                #Now copy GFs in, this will OVERWRITE EXISTING GFs of the same name
                flist=glob.glob(dirs[0]+'/*')
                for k in range(len(flist)):
                    copy(flist[k],outgreen)
                #Cleanup
                rmtree(dirs[0])
                gc.collect()
            else:  #Static GFs
                copy('staticgf',green_path+'static/'+model_name+'.static.'+strdepth+'.sub'+subfault)
                #Cleanup
                remove('staticgf')     
        #How long was I working for?
        toc=time.time()
        print('GFs computed in '+str((toc-tic)/60)+' minutes...')
    else:
        print('GFs not necessary when using an elastic halfspace, exiting make_green')


def make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
            hot_start,dk,pmin,pmax,kmax,ncpus,insar=False,okada=False):
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
    from numpy import arange,savetxt,genfromtxt
    from os import path,makedirs,environ
    from shlex import split
    import subprocess
    
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    #Load source model for station-event distance computations
    source=genfromtxt(fault_file)
    #Create all output folders
    for k in range(len(source)):
        strdepth='%.4f' % source[k,3]
        subfault=str(k+1).rjust(4,'0')
        if static==0 and tsunami==False:
            subfault_folder=green_path+'dynamic/'+model_name+'_'+strdepth+'.sub'+subfault
            if path.exists(subfault_folder)==False:
                #It doesn't, make it, don't be lazy
                makedirs(subfault_folder)               
#        elif static==0 and tsunami==True: #Tsunami GFs
#            subfault_folder=green_path+'tsunami/'+model_name+'_'+strdepth+'.sub'+subfault
#            if path.exists(subfault_folder)==False:
#                #It doesn't, make it, don't be lazy
#                makedirs(subfault_folder)
        elif static==1 and tsunami==True: #Tsunami GFs
            subfault_folder=green_path+'tsunami/'+model_name+'_'+strdepth+'.sub'+subfault
            if path.exists(subfault_folder)==False:
                #It doesn't, make it, don't be lazy
                makedirs(subfault_folder)
    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%.8f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    #Make mpi system call
    print("MPI: Starting GFs computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'
    if static==1 and okada==True:
        print('Static Okada solution requested, no need to run GFs...')
        pass
    else:
        mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_green '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(dt)+' '+str(NFFT)+' '+str(static)+' '+str(dk)+' '+str(pmin)+' '+str(pmax)+' '+str(kmax)+' '+str(tsunami)+' '+str(insar)
        print(mpi)
        mpi=split(mpi)
        p=subprocess.Popen(mpi)
        p.communicate()
        
        
        
        
        
def make_parallel_teleseismics_green(home,project_name,station_file,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus,hot_start=0):
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
    from numpy import arange,savetxt,genfromtxt
    from os import path,makedirs,environ
    from shlex import split
    import subprocess
    
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    #Load source model for station-event distance computations
    source=genfromtxt(fault_file)
   
    #Create all output folders
    for k in range(len(source)):
        
        strdepth='%.4f' % source[k,3]
        subfault=str(k+1).rjust(4,'0')
        
        subfault_folder=green_path+'dynamic/'+model_name+'_'+strdepth+'.sub'+subfault
        if path.exists(subfault_folder)==False:
            #It doesn't, make it, don't be lazy
            makedirs(subfault_folder)               

    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%.8f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    
    #Make mpi system call
    print("MPI: Starting GFs computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'


    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_teleseismics_green '+home+' '+project_name+' '+str(time_epi)+' '+station_file+' '+model_name+' '+teleseismic_vel_mod+' '+str(endtime)
    print(mpi)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
        
        


   

#Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,tsunami,beta,
                    hot_start,time_epi,impulse=False,okada=False,mu=45e9,insar=False):
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

        
    OUT:
        Nothing
    '''
    from mudpy import green
    import datetime
    from numpy import loadtxt
    import gc
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    logpath=home+project_name+'/logs/'
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Now compute synthetics please, one sub fault at a time
    for k in range(hot_start,source.shape[0]):
        print('ksource = ' + str(k))
        subfault=str(k+1).rjust(4,'0')
        log=green.run_syn(home,project_name,source[k,:],station_file,green_path,model_name,integrate,static,tsunami,
                subfault,time_epi,beta,impulse,okada,mu,insar=insar)
        f=open(logpath+'make_synth.'+now+'.log','a')
        f.write(log)
        f.close()
        gc.collect()
        
        
#Now make synthetics for source/station pairs
def make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,
                             tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse=False,insar=False,okada=False,mu=45e9):
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
    fault_file=home+project_name+'/data/model_info/'+fault_name
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Create individual source files
    for k in range(ncpus):
        i=arange(k+hot_start,len(source),ncpus)
        mpi_source=source[i,:]
        fmt='%d\t%10.6f\t%10.6f\t%.8f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f'
        savetxt(home+project_name+'/data/model_info/mpi_source.'+str(k)+'.fault',mpi_source,fmt=fmt)
    #Make mpi system call
    print("MPI: Starting synthetics computation on", ncpus, "CPUs\n")
    mud_source=environ['MUD']+'/src/python/mudpy/'
    mpi='mpiexec -n '+str(ncpus)+' python '+mud_source+'parallel.py run_parallel_synthetics '+home+' '+project_name+' '+station_file+' '+model_name+' '+str(integrate)+' '+str(static)+' '+str(quasistatic2dynamic)+' '+str(tsunami)+' '+str(time_epi)+' '+str(beta)+' '+str(custom_stf)+' '+str(impulse)+' '+str(insar)+' '+str(okada)+' '+str(mu)+' '+str(NFFT)+' '+str(dt)
    print(mpi)
    mpi=split(mpi)
    p=subprocess.Popen(mpi)
    p.communicate()
        
       
        
         
#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,tgf_file,fault_name,model_name,
        dt,tsun_dt,NFFT,tsunNFFT,green_flag,synth_flag,dk,pmin,
        pmax,kmax,beta,time_epi,hot_start,ncpus,custom_stf,quasistatic2dynamic=0,
        impulse=False):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,where,loadtxt,shape,floor
    from os import remove
    from gc import collect
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,dtype='U')
    GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],dtype='f8')
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    source=genfromtxt(fault_file)
    num_faults=shape(source)[0]
    if num_faults/ncpus < 2:
        ncpus=int(floor(num_faults/2.))
        print('Cutting back to ' + str(ncpus) + ' cpus for ' + str(num_faults) + ' subfaults')
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
            static=1
            tsunami=False
            insar=False
            make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus,insar)
        i=where(GF[:,3]==1)[0]
        if len(i)>0 : #displ waveform
            print('Displacememnt GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=str(stations[i[k]])+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=0
            tsunami=False
            make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        i=where(GF[:,4]==1)[0]
        if len(i)>0 : #vel waveform
            print('Velocity GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            static=0
            tsunami=False
            make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        if tgf_file!=None: #Tsunami
            print('Seafloor displacement GFs requested...')
#            static=0
            static=1
            tsunami=True
            station_file=tgf_file
            make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus)
        i=where(GF[:,6]==1)[0]
        if len(i)>0: #InSAR LOS
            print('InSAR GFs requested...')
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)): #Write temp .sta file
                #out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                out='%s\t%.8f\t%.8f\n' % (stations[i[k]],GF[i[k],0],GF[i[k],1])
                f.write(out)
            f.close()
            static=1
            tsunami=False
            insar=True
            make_parallel_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,tsunami,
                        hot_start,dk,pmin,pmax,kmax,ncpus,insar)
            collect()   
    #Synthetics are computed  one station at a time
    if synth_flag==1:
        #Paralell processing
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
            integrate=0
            static=1
            tsunami=False
            insar=False
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse,insar)
        #Decide which synthetics are required
        i=where(GF[:,3]==1)[0]
        if len(i)>0: #dispalcement waveform
            print('Displacement synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=str(stations[i[k]])+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=1
            static=0
            if tgf_file==None: # I am computing for stations on land
                tsunami=False
            else: #I am computing seafloor deformation for tsunami GFs, eventaully
                tsunami=True
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse)
        #Decide which synthetics are required
        i=where(GF[:,4]==1)[0]
        if len(i)>0: #velocity waveform
            print('Velocity synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=0
            static=0
            if tgf_file==None: # I am computing for stations on land
                tsunami=False
            else: #I am computing seafloor deformation for tsunami GFs, eventaully
                tsunami=True
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse)
        #Decide which synthetics are required
        i=where(GF[:,5]==1)[0]
        if len(i)>0: #tsunami waveform
            print('Tsunami synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            integrate=1
            static=1
            tsunami=True
            station_file=tgf_file
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse)
        #Decide which synthetics are required
        i=where(GF[:,6]==1)[0]
        if len(i)>0: # InSAR LOS
            print('InSAR synthetics requested')
            #Make dummy station file
            f=open(home+project_name+'/data/station_info/'+station_file,'w')
            for k in range(len(i)):
                #out=stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                out='%s\t%.8f\t%.8f\n' % (stations[i[k]],GF[i[k],0],GF[i[k],1])
                f.write(out)
            f.close()
            integrate=0
            static=1
            tsunami=False
            insar=True
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse,insar)
    
                   

    if quasistatic2dynamic == 1: #Conver static offsets to dyanmic step functions
    
        i = where(GF[:,2] == 1)[0]  #read from statics column in gflist
        if len(i) > 0: 
            
            print('quasistatic2dynamic requested')
            #Make dummy station file
            f = open(home+project_name+'/data/station_info/'+station_file,'w')
            
            for k in range(len(i)):
                out = stations[i[k]]+'\t'+repr(GF[i[k],0])+'\t'+repr(GF[i[k],1])+'\n'
                f.write(out)
            f.close()
            
            integrate=0
            static=0
            tsunami=False
            insar=True
            make_parallel_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,quasistatic2dynamic,tsunami,beta,hot_start,time_epi,ncpus,custom_stf,NFFT,dt,impulse,insar)
            
        else:
            
            print('ERROR: No sites in gflist file')
        



#Compute GFs for the ivenrse problem            
def teleseismicGFs(home,project_name,GF_list_teleseismic,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt,shape,floor
    from os import remove
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list_teleseismic
    stations=genfromtxt(gf_file,usecols=0,dtype='U')
    lonlat=genfromtxt(gf_file,usecols=[1,2,])
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    source=genfromtxt(fault_file)
    num_faults=shape(source)[0]
    
    if num_faults/ncpus < 2:
        ncpus=int(floor(num_faults/2.))
        print('Cutting back to ' + str(ncpus) + ' cpus for ' + str(num_faults) + ' subfaults')
    # GFs can be computed all at the same time
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    
    print('Teleseismic GFs requested...')
    f=open(home+project_name+'/data/station_info/'+station_file,'w')
    for k in range(len(stations)): #Write temp .sta file
        out=stations[k]+'\t'+repr(lonlat[k,0])+'\t'+repr(lonlat[k,1])+'\n'
        f.write(out)
    f.close()

    make_parallel_teleseismics_green(home,project_name,station_file,fault_name,model_name,teleseismic_vel_mod,time_epi,endtime,ncpus)
 




                              
                                                        
def run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,reg_spatial,reg_temporal,nfaults,beta,decimate,bandpass,
                solver,bounds,weight=False,Ltype=2,target_moment=None,data_vector=None,weights_file=None,
                onset_file=None,GOD_inversion=False):
    '''
    Assemble G and d, determine smoothing and run the inversion
    '''
    from mudpy import inverse as inv
    from mudpy.forward import get_mu_and_area
    from numpy import zeros,dot,array,squeeze,expand_dims,empty,tile,eye,ones,arange,load,size,genfromtxt
    from numpy import where,sort,r_,diag
    from numpy.linalg import lstsq
    from scipy.linalg import norm
    from scipy.sparse import csr_matrix as sparse
    from scipy.optimize import nnls
    from datetime import datetime
    import gc
    from matplotlib import path
    
    

    t1=datetime.now()
    #Get data vector
    if data_vector==None:
        d=inv.getdata(home,project_name,GF_list,decimate,bandpass=bandpass)
    else:
        d=load(data_vector) 
    #Get GFs
    G=inv.getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,decimate,bandpass,onset_file=onset_file)
    
    
    
    # Force faults inside a polygon to be zero (not contribute to inversion, 
    # useful for testing sensitivites)
    # print(' DANGER WILL ROBINSON: Forcing faults in polygon to have GFs = 0')
    
    # print('Keep Zone 1, 2 and 3')
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone1.txt')
    # zone_poly1 = path.Path(p)
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone2.txt')
    # zone_poly2 = path.Path(p)
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone3.txt')
    # zone_poly3 = path.Path(p)
    
    # # #Find indices of faults INSIDE polygon
    # fault_geometry = genfromtxt('/Users/dmelgarm/Slip_inv/Coalcoman2022/output/inverse_models/_old/models/gnss_sm_dart_insar_v3.0.0005.inv')
    # # The above is horrible form so I should explain. It's easier to read an alreayd existing ivnersion that has 
    # # the N windows already assigned and find the faults that are in the polygon that way
    # # I know. It's gross. Don't look at me that way.
    # izone = where((zone_poly1.contains_points(fault_geometry[:,1:3])== False) &
    #               (zone_poly2.contains_points(fault_geometry[:,1:3])== False) &  
    #               (zone_poly3.contains_points(fault_geometry[:,1:3])== False))[0]

    # # #Double indices because of ss and ds coordiante system
    # iss_zone = 2*izone
    # ids_zone = 2*izone + 1

    # #Zero out those GFs
    # G[:,iss_zone] = 0
    # G[:,ids_zone] = 0
    
 
    # print('Keep Zone 2 and 3')
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone2.txt')
    # zone_poly2 = path.Path(p)
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone3.txt')
    # zone_poly3 = path.Path(p)
    
    # # #Find indices of faults INSIDE polygon
    # fault_geometry = genfromtxt('/Users/dmelgarm/Slip_inv/Coalcoman2022/output/inverse_models/_old/models/gnss_sm_dart_insar_v3.0.0005.inv')
    # # The above is horrible form so I should explain. It's easier to read an alreayd existing ivnersion that has 
    # # the N windows already assigned and find the faults that are in the polygon that way
    # # I know. It's gross. Don't look at me that way.
    # izone = where((zone_poly2.contains_points(fault_geometry[:,1:3])== False) &  
    #               (zone_poly3.contains_points(fault_geometry[:,1:3])== False))[0]

    # # #Double indices because of ss and ds coordiante system
    # iss_zone = 2*izone
    # ids_zone = 2*izone + 1

    # #Zero out those GFs
    # G[:,iss_zone] = 0
    # G[:,ids_zone] = 0
    
    # print('Keep Zone 3')
    # p = genfromtxt('/Users/dmelgarm/Coalcoman2022/kml/zone3.txt')
    # zone_poly = path.Path(p)
    
    # #Find indices of faults INSIDE polygon
    # fault_geometry = genfromtxt('/Users/dmelgarm/Slip_inv/Coalcoman2022/output/inverse_models/_old/models/gnss_sm_dart_insar_v3.0.0005.inv')
    # # The above is horrible form so I should explain. It's easier to read an alreayd existing 
    # # ivnersion that has the N windows already assigned and find the faults that are in the 
    # # polygon that way. I know. It's gross. Don;t look at me that way.
    # #izone = where(zone_poly.contains_points(fault_geometry[:,1:3])==True)[0]  # <- if IN poly
    # izone = where(zone_poly.contains_points(fault_geometry[:,1:3])==False)[0]  # <- if OUT of poly

    # #Double indices because of ss and ds coordiante system
    # iss_zone3 = 2*izone
    # ids_zone3 = 2*izone + 1
    
    # #Zero out those GFs
    # G[:,iss_zone3] = 0
    # G[:,ids_zone3] = 0
    
    
    
    
    #print(' DANGER WILL ROBINSON: Forcing faults in polygon to have GFs = 0')
    

    
    #Find indices of faults INSIDE polygon
    # fault_geometry = genfromtxt('/Users/dmelgarm/Slip_inv/Turkey2023_M7.5_Surgu2/output/inverse_models/models/gnss_vr4.0.0010.inv')
    # # The above is horrible form so I should explain. It's easier to read an alreayd existing ivnersion that has 
    # # the N windows already assigned and find the faults that are in the polygon that way
    # # I know. It's gross. Don't look at me that way.

        
    # # izone = where((fault_geometry[:,3]<0.3) | (fault_geometry[:,1]>37.8))[0]
    # izone = where((fault_geometry[:,3]<0.5))[0]

    # #Double indices because of ss and ds coordiante system
    # iss_zone = 2*izone
    # ids_zone = 2*izone + 1

    # #Zero out those GFs
    # G[:,iss_zone] = 0
    # G[:,ids_zone] = 0
    

    
    #######    END POLY FILTER STUFF
    
    
    
    
    print(G.shape)                
    gc.collect()
    #Get data weights
    if weight==True:
        print('Applying data weights')
        if weights_file==None:
            w=inv.get_data_weights(home,project_name,GF_list,d,decimate)
        else:  # Remember weights are "uncertainties" alrger value is less trustworthy
            w = genfromtxt(weights_file)
            w = 1/w
        
        #apply data weights
        wd=w*d.squeeze()
        
        #get norm after applying weights and normalize weghted data vector
        data_norm = norm(wd)
        wd = wd/data_norm
        
        #reshape for inversion
        wd=expand_dims(wd,axis=1)
        
        #Apply weights to left hand side of the equation (the GFs)
        W=empty(G.shape)
        #W=tile(w,(G.shape[1],1)).T #why this and not a diagonal matrix??? ANSWER: they have the same effect, don;t rememebr why I chose to do it this way
        W=diag(w)
        
        #Normalization effect
        W = W/data_norm
        WG=W.dot(G)

        #Clear up extraneous variables
        W=None
        w=None
        #Define inversion quantities
        x=WG.transpose().dot(wd)
        print('Computing G\'G')
        K=(WG.T).dot(WG)
    else:
        #Define inversion quantities if no weighted
        x=G.transpose().dot(d)
        print('Computing G\'G')
        K=(G.T).dot(G)
    #Get regularization matrices (set to 0 matrix if not needed)
    static=False #Is it jsut a static inversion?
    if size(reg_spatial)>1:
        if Ltype==2: #Laplacian smoothing
            Ls=inv.getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
        elif Ltype==0: #Tikhonov smoothing
            N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
            Ls=eye(N) 
        elif Ltype==3:  #moment regularization
            N=nfaults[0]*nfaults[1]*num_windows*2 #Get total no. of model parameters
            Ls=ones((1,N))
            #Add rigidity and subfault area
            mu,area=get_mu_and_area(home,project_name,fault_name,model_name)
            istrike=arange(0,N,2)
            Ls[0,istrike]=area*mu
            idip=arange(1,N,2)
            Ls[0,idip]=area*mu
            #Modify inversion quantities
            x=x+Ls.T.dot(target_moment)
        else:
            print('ERROR: Unrecognized regularization type requested')
            return
        Ninversion=len(reg_spatial)
    else:
        Ls=zeros(K.shape)
        reg_spatial=array([0.])
        Ninversion=1
    if size(reg_temporal)>1:
        Lt=inv.getLt(home,project_name,fault_name,num_windows)
        Ninversion=len(reg_temporal)*Ninversion
    else:
        Lt=zeros(K.shape)
        reg_temporal=array([0.])
        static=True
    #Make L's sparse
    Ls=sparse(Ls)
    Lt=sparse(Lt)
    #Get regularization tranposes for ABIC
    LsLs=Ls.transpose().dot(Ls)
    LtLt=Lt.transpose().dot(Lt)
    #inflate
    Ls=Ls.todense()
    Lt=Lt.todense()
    LsLs=LsLs.todense()
    LtLt=LtLt.todense()
    #off we go
    dt=datetime.now()-t1
    print('Preprocessing wall time was '+str(dt))
    print('\n--- RUNNING INVERSIONS ---\n')
    ttotal=datetime.now()
    kout=0
    for kt in range(len(reg_temporal)):
        for ks in range(len(reg_spatial)):
            t1=datetime.now()
            lambda_spatial=reg_spatial[ks]
            lambda_temporal=reg_temporal[kt]
            print('Running inversion '+str(kout+1)+' of '+str(Ninversion)+' at regularization levels: ls ='+repr(lambda_spatial)+' , lt = '+repr(lambda_temporal))
            if static==True: #Only statics inversion no Lt matrix
                Kinv=K+(lambda_spatial**2)*LsLs
                Lt=eye(len(K))
                LtLt=Lt.T.dot(Lt)
            else: #Mixed inversion
                Kinv=K+(lambda_spatial**2)*LsLs+(lambda_temporal**2)*LtLt
            if solver.lower()=='lstsq':
                sol,res,rank,s=lstsq(Kinv,x)
            elif solver.lower()=='nnls':
                x=squeeze(x.T)
                try:
                    sol,res=nnls(Kinv,x)
                except:
                    print('+++ WARNING: No solution found, writting zeros.')
                    sol=zeros(G.shape[1])
                x=expand_dims(x,axis=1)
                sol=expand_dims(sol,axis=1)
            else:
                print('ERROR: Unrecognized solver \''+solver+'\'')

            #Force faults outside a polygon to be zero
#            print('WARNING: Using fault polygon to force solutions to zero')
#            #load faulkt
#            fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
#            polygon=genfromtxt('/Users/dmelgarm/Oaxaca2020/etc/zero_fault.txt')
#            polygon=path.Path(polygon)
#            i=where(polygon.contains_points(fault[:,1:3])==False)[0]
#            i=sort(r_[i*2,i*2+1])
#            N=nfaults[0]*2
#            i=r_[i,i+N,i+2*N,i+3*N]
#            sol[i]=0
            
            #Compute synthetics
            ds=dot(G,sol)
            
            #Get stats
            L2,Lmodel=inv.get_stats(Kinv,sol,x,Ls)
            VR,L2data=inv.get_VR(home,project_name,GF_list,sol,d,ds,decimate,WG,wd)
            #VR=inv.get_VR(WG,sol,wd)
            #ABIC=inv.get_ABIC(WG,K,sol,wd,lambda_spatial,lambda_temporal,Ls,LsLs,Lt,LtLt)
            ABIC=inv.get_ABIC(G,K,sol,d,lambda_spatial,lambda_temporal,Ls,LsLs,Lt,LtLt)
            #Get moment
            Mo,Mw=inv.get_moment(home,project_name,fault_name,model_name,sol)
            #If a rotational offset was applied then reverse it for output to file
            if beta !=0:
                sol=inv.rot2ds(sol,beta)
            #Write log
            inv.write_log(home,project_name,run_name,kout,rupture_speed,num_windows,lambda_spatial,lambda_temporal,beta,
                L2,Lmodel,VR,ABIC,Mo,Mw,model_name,fault_name,G_name,GF_list,solver,L2data)
            #Write output to file
            if GOD_inversion==True:
                num=str(kout).rjust(4,'0')
                np.save(home+project_name+'/output/inverse_models/'+run_name+'.'+num+'.syn.npy',ds)
                inv.write_synthetics_GOD(home,project_name,run_name,GF_list,ds,kout,decimate)
            else:
                inv.write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,kout,decimate)
            inv.write_model(home,project_name,run_name,fault_name,model_name,rupture_speed,num_windows,epicenter,sol,kout,onset_file=onset_file)
            kout+=1
            dt1=datetime.now()-t1
            dt2=datetime.now()-ttotal
            print('... inversion wall time was '+str(dt1)+', total wall time elapsed is '+str(dt2))
