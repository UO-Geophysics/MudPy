"""
Created on Wed Dec 19 14:28:43 2018

@author: dmelgarm
"""

def get_station_delays(station_coords,sources,velmod='PREM',phase_list=['s','S']):
    '''
    Given an array of station corodiantes and sources calculate the travel time
    from each source toe ach site.
    
    velmod is the FULL path to the .npz file used by TauPy
    '''
    
    from obspy.taup import TauPyModel
    from numpy import zeros
    from obspy.geodetics.base import locations2degrees

    model = TauPyModel(model=velmod)
    
    #Initalize output variabe
    delay_time=zeros((len(sources),len(station_coords)))
    
    #Loop over sources
    for ksource in range(len(sources)):
        
        print('Working on source %d of %d ' % (ksource,len(sources)))
        
        #loop over sites
        for ksite in range(len(station_coords)):
        
            distance_in_degrees=locations2degrees(station_coords[ksite,1],station_coords[ksite,0],
                                                      sources[ksource,2],sources[ksource,1])
            
            arrivals = model.get_travel_times(source_depth_in_km=sources[ksource,3],
                                  distance_in_degree=distance_in_degrees,phase_list=phase_list)
            
            delay_time[ksource,ksite]=arrivals[0].time
            
    return delay_time






def run_travel_times(home,project_name,model_name,fault_name,stations_list,
                     traveltimes_name):
    '''
    This function runs everything to create the trasvel tiems table and 
    save it as an npz file
    '''
    
    from numpy import genfromtxt,savez
    
    station_coords=genfromtxt(home+project_name+'/data/station_info/'+stations_list,usecols=[1,2])
    velmod=model_name
    sources=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    
    #Call travel times function
    travel_times=get_station_delays(station_coords,sources,velmod=velmod,phase_list=['s','S'])
    
    #save to npz file
    savez(home+project_name+'/travel_times/'+traveltimes_name+'.npz',travel_times)
    
    
def backproject(home,project_name,run_name,model_name,fault_name,stations_list,
                     traveltimes_name,time_epi,Tmax,stack_order=4):
    '''
    Form the waveform stacks and back rooject toe ach source as a function of
    time. The output will be a stream object with one trace that contains the 
    stack at each source point
    '''
    
    from numpy import genfromtxt,load
    from obspy import read,Stream,Trace
    from mudpy.strong_motion import Nstack
    

    
    #laod npz file with travel times
    t=load(home+project_name+'/travel_times/'+traveltimes_name+'.npz')
    delay_times=t[t.files[0]]
    
    #load sources
    sources=genfromtxt(home+project_name+'/data/model_info/'+fault_name)

    #get station paths
    station_paths=genfromtxt(home+project_name+'/data/station_info/'+stations_list,usecols=3,dtype='S')
    
    #read data
    for ksite in range(len(station_paths)):
        
        if ksite==0:
            st=read(station_paths[0])
        else:
            st+=read(station_paths[ksite])
            
        #trim and remove baseline
        st[ksite].trim(starttime=time_epi,endtime=time_epi+Tmax,pad=True,fill_value=0)
            
    #Form stack for each source point
    stack=Stream(Trace())
    
    for ksource in range(len(sources)):
        
        #This gets reinitalized for each new source
        st_for_stack=st.copy()
        
        print('... working on stack %d of %d ' %(ksource,len(sources)))
        
        for ksite in range(len(st)):
            
            #get current delay time
            dt=delay_times[ksource,ksite]
            
            #delay waveform by cropping the right amount of samples
            N_crop_samples=int(dt/st[ksite].stats.delta)
            data=st[ksite].data
            data=data[N_crop_samples:-1]
            
            #place back in a trace object
            st_for_stack[ksite].data=data
            
            #Trim to pad ends with zeros
            st_for_stack[ksite].trim(endtime=st[0].stats.endtime,pad=True,fill_value=0)
            
        #Form the stack
        if ksource==0:
            stack[0]=Nstack(st_for_stack,N=stack_order,normalize=True)
        else:
            stack+=Nstack(st_for_stack,N=stack_order,normalize=True)
            
            
    stack.write(home+project_name+'/output/models/'+run_name+'.stacks.mseed',fomrat='MSEED')
            
            
        
            
def make_frames(home,project_name,run_name,fault_name,dt_frame=1.0):
    '''
    Read the mseed file with the stacks and make an array witht he snapshots
    every dt_frames seconds
    '''
    
    from numpy import genfromtxt,zeros,arange,where
    from obspy import read
    
    
    #read stacks
    stacks=read(home+project_name+'/output/models/'+run_name+'.stacks.mseed')
    
    #read sources
    sources=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    
    #How many frames?
    total_time=stacks[0].stats.endtime-stacks[0].stats.starttime
    time_frames=arange(0,total_time+1e-6,dt_frame)
    
    #make frames
    frames=zeros((len(sources),len(time_frames)))
    
    for ktime in range(len(time_frames)):
    
        print('Working on time slice %d of %d' % (ktime,len(time_frames)))
        
        #find correct time
        i=where(stacks[0].times()==time_frames[ktime])[0]
        
        for ksource in range(len(sources)):
            
                frames[ksource,ktime]=stacks[ksource].data[i]
            
    