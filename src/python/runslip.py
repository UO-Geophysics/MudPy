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
        makedirs(proj_dir+'data/statics')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'forward_models')
        makedirs(proj_dir+'output/inverse_models')
        makedirs(proj_dir+'output/inverse_models/statics')
        makedirs(proj_dir+'output/inverse_models/waveforms')
        makedirs(proj_dir+'output/inverse_models/models')
        makedirs(proj_dir+'output/forward_models')
        makedirs(proj_dir+'logs')


#Extract fault geometry from rupture file
def rupt2fault(home,project_name,rupture_name):
    '''
    Make fault file from user provided forward model rupture file
    '''
    from numpy import loadtxt,savetxt,c_
    
    print 'Assembling fault file from rupture file'
    rupt=loadtxt(home+project_name+'/forward_models/'+rupture_name,ndmin=2)
    fault=c_[rupt[:,0],rupt[:,1],rupt[:,2],rupt[:,3],rupt[:,4],rupt[:,5],rupt[:,6],rupt[:,7],rupt[:,10],rupt[:,11]]
    savetxt(home+project_name+'/data/model_info/'+rupture_name.split('.')[0]+'.fault', \
            fault,fmt='%i\t%.6f\t%.6f\t%.6f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f')





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
    for k in range(source.shape[0]):
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
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,beta,
                    hot_start,coord_type,time_epi):
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
        log=green.run_syn(home,project_name,source[k,:],station_file,green_path,model_name,integrate,static,
                subfault,coord_type,time_epi,beta)
        f=open(logpath+'make_synth.'+now+'.log','a')
        f.write(log)
        f.close()
       
        
         
#Compute GFs for the ivenrse problem            
def inversionGFs(home,project_name,GF_list,fault_name,model_name,dt,NFFT,coord_type,
                green_flag,synth_flag,dk,pmin,pmax,kmax,beta,time_epi,hot_start):
    '''
    This routine will read a .gflist file and compute the required GF type for each station
    '''
    from numpy import genfromtxt
    from os import remove
    from gc import collect
    
    #Read in GFlist and decide what to compute
    gf_file=home+project_name+'/data/station_info/'+GF_list
    stations=genfromtxt(gf_file,usecols=0,skip_header=1,dtype='S6')
    GF=genfromtxt(gf_file,usecols=[1,2,3,4,5,6,7],skip_header=1,dtype='f8')
    #Now do one station at a time
    station_file='temp.sta'
    try:
        remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
    except:
        pass
    for k in range(hot_start,len(stations)):
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
            collect()
        if synth_flag==1:
            #Decide which synthetics are required
            if GF[k,2]==1: #Static offset
                integrate=0
                static=1
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,beta,hot_start,coord_type,time_epi)
            if GF[k,3]==1: #dispalcement waveform
                integrate=1
                static=0
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,beta,hot_start,coord_type,time_epi)
            if GF[k,4]==1: #velocity waveform
                integrate=0
                static=0
                make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,beta,hot_start,coord_type,time_epi)
            if GF[k,5]==1: #tsunami waveform
                pass
            if GF[k,6]==1: #strain offsets
                pass
            
    remove(home+project_name+'/data/station_info/'+station_file) #Cleanup
                    
                                
                                                        
def run_inversion(home,project_name,run_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,coord_type,reg_spatial,reg_temporal,nfaults,beta,decimate,lowpass,
                solver,bounds):
    '''
    Assemble G and d, determine smoothing and run the inversion
    '''
    import inverse as inv
    from numpy import zeros,dot,array,squeeze,expand_dims,empty,tile
    from numpy.linalg import lstsq
    from scipy.sparse import csr_matrix as sparse
    from scipy.optimize import nnls
    from datetime import datetime
    
    t1=datetime.now()
    #Get data vector
    d=inv.getdata(home,project_name,GF_list,decimate,lowpass)
    #Get GFs
    G=inv.getG(home,project_name,fault_name,model_name,GF_list,G_from_file,G_name,epicenter,
                rupture_speed,num_windows,coord_type,decimate)
    #Get data weights
    w=inv.get_data_weights(home,project_name,GF_list,d,decimate)
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
    print 'Computing G\'G'
    K=(WG.T).dot(WG)
    #Get regularization matrices (set to 0 matrix if not needed)
    if type(reg_spatial)!=bool:
        Ls=inv.getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
        Ninversion=len(reg_spatial)
    else:
        Ls=zeros(K.shape)
        reg_spatial=array([0.])
        Ninversion=1
    if type(reg_temporal)!=bool:
        Lt=inv.getLt(home,project_name,fault_name,num_windows)
        Ninversion=len(reg_temporal)*Ninversion
    else:
        Lt=zeros(K.shape)
        reg_temporal=array([0.])
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
    kout=0
    dt=datetime.now()-t1
    print 'Preprocessing wall time was '+str(dt)
    print '\n--- RUNNING INVERSIONS ---\n'
    ttotal=datetime.now()
    for kt in range(len(reg_temporal)):
        for ks in range(len(reg_spatial)):
            t1=datetime.now()
            lambda_spatial=reg_spatial[ks]
            lambda_temporal=reg_temporal[kt]
            print 'Running inversion '+str(kout+1)+' of '+str(Ninversion)+' at regularization levels: ls ='+repr(lambda_spatial)+' , lt = '+repr(lambda_temporal)
            Kinv=K+(lambda_spatial**2)*LsLs+(lambda_temporal**2)*LtLt
            if solver.lower()=='lstsq':
                sol,res,rank,s=lstsq(Kinv,x)
            elif solver.lower()=='nnls':
                x=squeeze(x.T)
                sol,res=nnls(Kinv,x)
                x=expand_dims(x,axis=1)
                sol=expand_dims(sol,axis=1)
            else:
                print 'ERROR: Unrecognized solver \''+solver+'\''
            #Compute synthetics
            ds=dot(G,sol)
            #Get stats
            L2,Lmodel=inv.get_stats(Kinv,sol,x)
            VR=inv.get_VR(G,sol,d)
            ABIC=inv.get_ABIC(WG,K,sol,wd,lambda_spatial,lambda_temporal,Ls,LsLs,Lt,LtLt)
            #Get moment
            Mo,Mw=inv.get_moment(home,project_name,fault_name,model_name,sol)
            #If a rotational offset was applied then reverse it for output to file
            if beta !=0:
                sol=inv.rot2ds(sol,beta)
            #Write log
            inv.write_log(home,project_name,run_name,kout,rupture_speed,num_windows,lambda_spatial,lambda_temporal,beta,
                L2,Lmodel,VR,ABIC,Mo,Mw,model_name,fault_name,G_name,GF_list,solver)
            #Write output to file
            inv.write_synthetics(home,project_name,run_name,GF_list,G,sol,ds,kout,decimate)
            inv.write_model(home,project_name,run_name,fault_name,model_name,rupture_speed,num_windows,epicenter,sol,kout)
            kout+=1
            dt1=datetime.now()-t1
            dt2=datetime.now()-ttotal
            print '... inversion wall time was '+str(dt1)+', total wall time elapsed is '+str(dt2)