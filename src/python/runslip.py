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
    from shutil import rmtree
    from os import path,makedirs
    clob='y'
    proj_dir=home+project_name+'/'
    if path.exists(proj_dir):  #Path exists, clobber?
        clob=raw_input('Project directory exists, clobber (y/n)?')
        if clob is'y' or clob is 'Y': #Clobber baby
            clob=raw_input('This will delete everything in this project directory, so, take a minute, think about it: clobber (y/n)?')
            if clob is 'y' or clob is 'Y':
                rmtree(proj_dir)
            else: #Leave direcory alone
                print 'Phew, almost shot yourself in the foot there didn\'t you?'
        else: #Leave direcory alone
            print 'Phew, almost shot yourself in the foot there didn\'t you?'
    if clob is 'y' or clob is 'Y':
        makedirs(proj_dir)
        #And make the subdirectories
        makedirs(proj_dir+'GFs')
        makedirs(proj_dir+'GFs/static')
        makedirs(proj_dir+'GFs/dynamic')
        makedirs(proj_dir+'GFs/matrices')
        makedirs(proj_dir+'data/waveforms')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'forward_models')
        makedirs(proj_dir+'output/inverse_models')
        makedirs(proj_dir+'output/forward_models')
        makedirs(proj_dir+'logs')


#Setup some book-keeping for the forward problem
def forward_setup(home,project_name,rupture_name):
    '''
    Make fault file from user provided forward model rupture file
    '''
    from numpy import loadtxt,savetxt,c_
    
    print '======= FORWARD MODELING ========'
    rupt=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    fault=c_[rupt[:,0],rupt[:,1],rupt[:,2],rupt[:,3],rupt[:,4],rupt[:,5],rupt[:,7],rupt[:,8]]
    savetxt(home+project_name+'/data/model_info/'+rupture_name.split('.')[0]+'.fault', \
            fault,fmt='%i\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.4f\t%.4f')





# Run green functions          
def make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start,coord_type,dk,pmin,pmax,kmax):
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
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        Nothing
    '''
    import time,glob,green
    from numpy import loadtxt
    from shutil import rmtree,copy
    from os import chdir,path,makedirs,remove
    from string import rjust
    import datetime
    
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
        log=green.run_green(source[k,:],station_file,model_name,dt,NFFT,static,coord_type,dk,pmin,pmax,kmax)
        #Write log
        f=open(logpath+'make_green.'+now+'.log','a')
        f.write(log)
        f.close()
        #Move to correct directory
        strdepth='%.4f' % source[k,3]
        subfault=rjust(str(k+1),4,'0')
        if static==0:
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
        else:  #Static GFs
            copy('staticgf',green_path+'static/'+model_name+'.static.'+strdepth+'.sub'+subfault)
            #Cleanup
            remove('staticgf')     
    #How long was I working for?
    toc=time.time()
    print 'GFs computed in '+str((toc-tic)/60)+' minutes...'




#Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start,coord_type,time_epi):
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
    import green,datetime
    from numpy import loadtxt
    from string import rjust
    
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
        subfault=rjust(str(k+1),4,'0')
        log=green.run_syn(source[k,:],station_file,green_path,model_name,integrate,static,subfault,coord_type,time_epi)
        f=open(logpath+'make_synth.'+now+'.log','a')
        f.write(log)
        f.close()
       
        
         
#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,fault_name,model_name,dt,NFFT,coord_type,green_flag,synth_flag,dk,pmin,pmax,kmax,time_epi):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt
    from os import remove
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='S6')
    GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
    #Now do one station at a time
    hot_start=0 #Used for the forward problem, doesn't do anything here
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    for k in range(len(stations)):
        #Make dummy station file
        out=stations[k]+'\t'+repr(GF[k,0])+'\t'+repr(GF[k,1])
        f=open(home+project_name+'/data/station_info/'+station_file,'w')
        f.write(out)
        f.close()
        if green_flag==1:
            #decide what GF computation is required for this station
            if GF[k,2]==1: #Static offset
                static=1
                make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start,coord_type,dk,pmin,pmax,kmax)
            if GF[k,3]==1 or GF[k,4]==1: #full waveform
                static=0
                make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT,static,hot_start,coord_type,dk,pmin,pmax,kmax)
            if GF[k,5]==1: #Tsunami (pending)
                pass
            if GF[k,6]==1: #strain (pending)
                pass
        if synth_flag==1:
            #Decide which synthetics are required
            if GF[k,2]==1: #Static offset
                integrate=0
                static=1
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start,coord_type)
            if GF[k,3]==1: #dispalcement waveform
                integrate=1
                static=0
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start,coord_type,time_epi)
            if GF[k,4]==1: #velocity waveform
                integrate=0
                static=0
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,hot_start,coord_type,time_epi)
            if GF[k,5]==1: #tsunami waveform
                pass
            if GF[k,6]==1: #strain offsets
                pass
    remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
                                
def run_inversion(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speeds,coord_type,bounds,regularization_type,regularization_parameter,nfaults):
    '''
    Assemble G and d, determine smoothing and run the inversion
    '''
    from inverse import getG,getdata,getL
    
    #Get data vector
    d=getdata(home,project_name,GF_list,rupture_speeds)
    #Get GFs
    G=getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speeds,coord_type)
    L=getL(home,project_name,fault_name,bounds,regularization_type,nfaults)
    s,W=get_data_weights(home,project_name,GF_list)

        
                


######                 Le tools undt le trinkets                         #######
                                

        
        
        
    
        