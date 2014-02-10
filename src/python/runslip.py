'''
Diego Melgar, 01/2014
Runtime file for forward modeling and inverse kinematic slip inversions
'''


#Initalize project folders
def init(home,project_name):
    '''
    Blerg
    '''
    from shutil import rmtree,copy
    from os import path,makedirs
    clob='y'
    proj_dir=home+project_name+'/'
    aux_dir=home+'aux/'
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
        makedirs(proj_dir+'data/waveforms')
        makedirs(proj_dir+'data/station_info')
        makedirs(proj_dir+'data/model_info')
        makedirs(proj_dir+'structure')
        makedirs(proj_dir+'plots')
        makedirs(proj_dir+'output/forward_models')
        makedirs(proj_dir+'output/inverse_models')
        makedirs(proj_dir+'output/synthetics')
        #Copy aux files to proper places
        copy(aux_dir+'init.fault',proj_dir+'data/model_info/') #Blank fault file
        copy(aux_dir+'init.sta',proj_dir+'data/station_info/') #Blank station file
        copy(aux_dir+'init.mod',proj_dir+'structure/') #Blank structure file

# Run green functions          
def make_green(home,project_name,station_file,fault_name,model_name,dt,NFFT):
    '''
    di Blergerson
    '''
    import time,glob,green
    from numpy import loadtxt
    from shutil import rmtree,copy
    from os import chdir,path,makedirs
    
    tic=time.time()
    model_path=home+project_name+'/structure/'
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file 
    fault_file=home+project_name+'/data/model_info/'+fault_name  
    chdir(model_path)
    #Load source model for station-event distance computations
    source=loadtxt(fault_file,ndmin=2)
    #Pass only unique source points into computation (avoid duplication for strike-slip and dip-slip components)
    #this is PENDING, use numpy.unique and numpy.intersect1d, computation is not that slow so I might not do this
    for k in range(source.shape[0]):
        green.run_green(source[k,:],station_file,model_name,dt,NFFT)
        strdepth='%.1f' % source[k,3]
        #Move results to GF dir
        dirs=glob.glob('*.mod_'+strdepth)
        #Where am I writting this junk too?
        outgreen=green_path+path.split(dirs[0])[1]
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
    #How long was I out?
    toc=time.time()
    print 'GFs computed in '+str((toc-tic)/60)+' minutes...'




#Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate):
    '''
    Blergarmon
    '''
    import green
    from numpy import loadtxt
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    green_dir=green_path+model_name
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Now synthetics please, one sub fault at a time
    for k in range(source.shape[0]):
        green.run_syn(source[k,:],station_file,green_dir,integrate)