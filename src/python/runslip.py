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
def make_green(home,project_name,model_name,min_depth,max_depth,delta_depth,min_distance,max_distance,delta_distance,dt,NFFT):
    '''
    di Blergerson
    '''
    import time,glob,green
    import numpy as np
    from shutil import rmtree,move
    from os import chdir
    
    if delta_depth!=0:
        source_depths=np.arange(min_depth,max_depth,delta_depth)
    else:
        source_depths=np.array([min_depth])
        delta_depth=1
        station_distances=np.arange(min_distance,max_distance,delta_distance) 
    tic=time.time()
    model_path=home+project_name+'/structure/'
    green_path=home+project_name+'/GFs/'   
    chdir(model_path)  
    for k in range(len(source_depths)):
        green.run_green(source_depths[k],station_distances,model_name,dt,NFFT)
    #Move to GF dir
    dirs=glob.glob('*.mod_*')
    for k in range(len(source_depths)):
        #Delete previous GFs and move in the new ones
        try:
            rmtree(green_path+dirs[k])
        except:
            pass
        move(dirs[k],green_path)
    toc=time.time()
    print 'GFs computed in '+str((toc-tic)/60)+' minutes...'

#Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,delta_depth,delta_distance,integrate):
    '''
    Blergarmon
    '''
    import green
    from numpy import loadtxt,round
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    green_dir=green_path+model_name
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Now synthetics please, one sub fault at a time
    for k in range(source.shape[0]):
        #Round source depth to nearest computed depth
        if delta_depth==0:
            delta_depth=1
        source[k,3]=round(source[k,3]/delta_depth)*delta_depth
        green.run_syn(source[k,:],station_file,delta_distance,green_dir,integrate)
    
    
    