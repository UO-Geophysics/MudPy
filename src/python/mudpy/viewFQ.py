'''
Module for plotting FakeQuakes related studd
'''
from matplotlib import rc
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

rc('font', **font)

from obspy.imaging.cm import pqlx

#Colormap
whitejet_dict = {'red': ((0., 1, 1),
                 (0.02, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
                'green': ((0., 1, 1),
                   (0.02, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
                'blue': ((0., 1, 1),
                  (0.02, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

pqlx_dict={'blue': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  1.  ,  1.  ),
            ( 0.35 ,  1.  ,  1.  ),
            ( 0.6 ,  1.  ,  1.  ),
            ( 0.7 ,  0.  ,  0.  ),
            ( 0.89 ,  0.  ,  0.  ),
            ( 1.  ,  0.  ,  0.  )), 
            'green': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  0.  ,  0.  ),
            ( 0.35 ,  0.  ,  0. ),
            ( 0.6 ,  1.  ,  1.  ),
            ( 0.7 ,  1.  ,  1.  ),
            ( 0.89 ,  1.  ,  1.  ),
            ( 1.  ,  0.  ,  0. )), 
            'red': (( 0.  ,  1.  ,  1.  ),
            ( 0.1,  1.  ,  1.  ),
            ( 0.22,  1.  ,  1.  ),
            ( 0.35 ,  0.  ,  0.  ),
            ( 0.6 ,  0.  ,  0.  ),
            ( 0.7 ,  0.  ,  0.  ),
            ( 0.89 ,  1.  ,  1.  ),
            ( 1.  ,  1.  ,  1.  ))}


rise_times_dict={u'blue': [(0.0, 1.0, 1.0),
                        (0.06, 1.0, 1.0),
                        (0.125, 0.5, 0.5),
                        (0.25, 0.1, 0.1),
                        (0.375, 0.0, 0.0),
                        (0.5, 0.15, 0.15),
                        (0.625, 0.5, 0.5),
                        (0.75, 0.75, 0.75),
                        (0.875, 0.5, 0.5),
                        (1.0, 0.0, 0.0)],
                        u'green': [(0.0, 1.0, 1.0),
                        (0.06, 1.0, 1.0),
                        (0.125, 0.9, 0.9),
                        (0.25, 0.75, 0.75),
                        (0.375, 0.5, 0.5),
                        (0.5, 0.25, 0.25),
                        (0.625, 0.2, 0.2),
                        (0.75, 0.15, 0.15),
                        (0.875, 0.15, 0.15),
                        (1.0, 0.0, 0.0)],
                        u'red': [(0.0, 1.0, 1.0),
                        (0.06, 1.0, 1.0),
                        (0.125, 0.9, 0.9),
                        (0.25, 0.9, 0.9),
                        (0.375, 0.9, 0.9),
                        (0.5, 1.0, 1.0),
                        (0.625, 0.6, 0.6),
                        (0.75, 0.3, 0.3),
                        (0.875, 0.15, 0.15),
                        (1.0, 0.0, 0.0)]}



def plot_KLslip(home,project_name,run_name,mesh_name,fudge=0.3):
    '''
    Make simple plot of slip model
    '''
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib import pyplot as plt
    from matplotlib import cm,colors
    from numpy import r_,genfromtxt,reshape,zeros,array
    from glob import glob
    from string import replace

    #Read mesh
    mesh=genfromtxt(home+project_name+'/data/model_info/'+mesh_name)
    
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

    rc('font', **font)
    
    whitejet = colors.LinearSegmentedColormap('whitejet',whitejet_dict,256)
    pqlx_colormap = colors.LinearSegmentedColormap('pqlx',pqlx_dict,256)
    rise_times_colormap = colors.LinearSegmentedColormap('rtimes',rise_times_dict,256)
    #Get faults I need to work on
    faults=glob(home+project_name+'/output/ruptures/*.rupt')
    logs=glob(home+project_name+'/output/ruptures/*.log')
    for kfault in range(len(faults)):
        print faults[kfault]
        #Read slip info
        fault=genfromtxt(faults[kfault])
        
        #Get info about fault
        f=open(logs[kfault],'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Hypocenter (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo_fault=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])   
            if 'Target magnitude' in line:
                target_Mw=float(line.split(':')[-1].split(' ')[-1])        
        
        #Plot
        fig, axarr = plt.subplots(1,3,figsize=(12,8))
        
        ax1=axarr[0]
        ax1.set_xlim([fault[:,1].min()-fudge,fault[:,1].max()+2])
        ax1.set_ylim([fault[:,2].min()-fudge,fault[:,2].max()+fudge])
    
        #Plot coasts (also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/coast.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax1.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.5,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
        #Plot cnational boundaries(also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/boundaries.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax1.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=1,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
    
    
        #Make patches
        patches = []
        for k in range(len(fault)):
            coordinates=reshape(r_[mesh[k,4:6],mesh[k,7:9],mesh[k,10:12],mesh[k,4:6]],(4,2))
            subfault = Polygon(coordinates, True)
            patches.append(subfault)
            
        p1 = PatchCollection(patches,cmap=pqlx_colormap,lw=0.1)
        colors = fault[:,9]
        p1.set_array(colors)
        ax1.add_collection(p1)
        labels = ax1.get_xticklabels() 
        for label in labels: 
            label.set_rotation(70) 
        #Plot stations, this is hard coded to Cascadia
        stations=genfromtxt(u'/Users/dmelgar/Cascadia/stations/cascadia_gps_small.txt',usecols=[1,2])
        ax1.scatter(stations[:,0],stations[:,1],marker='o',lw=0.8,c='#FF8C00',s=15)
        #hypocenter
        ax1.scatter(hypo_fault[0],hypo_fault[1],marker='*',lw=0.9,c='#9ACD32',s=300)
        ax1.set_title('Target Mw = %.2f\nActual Mw = %.2f' % (target_Mw,Mw))
        
        
        #### Plot rise times
        
        ax2=axarr[1]
        ax2.set_xlim([fault[:,1].min()-fudge,fault[:,1].max()+2])
        ax2.set_ylim([fault[:,2].min()-fudge,fault[:,2].max()+fudge])
    
        #Plot coasts (also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/coast.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax2.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.5,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
        #Plot cnational boundaries(also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/boundaries.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax2.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=1,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
    
    
        #Make patches
        patches = []
        for k in range(len(fault)):
            coordinates=reshape(r_[mesh[k,4:6],mesh[k,7:9],mesh[k,10:12],mesh[k,4:6]],(4,2))
            subfault = Polygon(coordinates, True)
            patches.append(subfault)
            
        p2 = PatchCollection(patches,cmap=rise_times_colormap,lw=0.1)
        colors = fault[:,7]
        p2.set_array(colors)
        ax2.add_collection(p2)
        labels = ax2.get_xticklabels() 
        for label in labels: 
            label.set_rotation(70) 
        #Plot stations, this is hard coded to Cascadia
        stations=genfromtxt(u'/Users/dmelgar/Cascadia/stations/cascadia_gps_small.txt',usecols=[1,2])
        ax2.scatter(stations[:,0],stations[:,1],marker='o',lw=0.8,c='#FF8C00',s=15)
        #hypocenter
        ax2.scatter(hypo_fault[0],hypo_fault[1],marker='*',lw=0.9,c='#9ACD32',s=300)
        #plt.title('Target Mw = %.2f\nActual Mw = %.2f' % (target_Mw,Mw))
        
        
        #### Plot Rupture onset times
        
        ax3=axarr[2]
        ax3.set_xlim([fault[:,1].min()-fudge,fault[:,1].max()+2])
        ax3.set_ylim([fault[:,2].min()-fudge,fault[:,2].max()+fudge])
    
        #Plot coasts (also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/coast.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax3.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.5,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
        #Plot cnational boundaries(also hard coded to cascadia for now)
        fcoast=open('/Users/dmelgar/Cascadia/coast/boundaries.txt','r')
        line=fcoast.readline()
        parsing_coast=True
        while parsing_coast==True:
            if '>' in line:
                complete_polygon=False
                count=0
                while complete_polygon==False:
                    line=fcoast.readline()
                    if '>' in line or line=='':
                        complete_polygon=True
                        ax3.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=1,c='k')
                    else:
                        new_coordinates=zeros((1,2))
                        new_coordinates[0,0]=float(line.split()[0])
                        new_coordinates[0,1]=float(line.split()[1])
                        if count==0:
                            coast_coordinates=new_coordinates.copy()
                            count+=1
                        else:
                            coast_coordinates=r_[coast_coordinates,new_coordinates]
                    if line=='':
                        parsing_coast=False
        #Done plotting coast
    
    
        #Make patches
        patches = []
        for k in range(len(fault)):
            coordinates=reshape(r_[mesh[k,4:6],mesh[k,7:9],mesh[k,10:12],mesh[k,4:6]],(4,2))
            subfault = Polygon(coordinates, True)
            patches.append(subfault)
            
        p3 = PatchCollection(patches,cmap=whitejet,lw=0.1)
        colors=fault[:,12]
        p3.set_array(colors)
        ax3.add_collection(p3)
        labels = ax3.get_xticklabels() 
        for label in labels: 
            label.set_rotation(70) 
        #Plot stations, this is hard coded to Cascadia
        stations=genfromtxt(u'/Users/dmelgar/Cascadia/stations/cascadia_gps_small.txt',usecols=[1,2])
        ax3.scatter(stations[:,0],stations[:,1],marker='o',lw=0.8,c='#FF8C00',s=15)
        #hypocenter
        ax3.scatter(hypo_fault[0],hypo_fault[1],marker='*',lw=0.9,c='#9ACD32',s=300)
        #plt.title('Target Mw = %.2f\nActual Mw = %.2f' % (target_Mw,Mw))
        
                
        #Colorbars
        cbar_ax1 = fig.add_axes([0.32, 0.2, 0.015, 0.6])
        plt.colorbar(p1, cax=cbar_ax1, label="Slip (m)")
        cbar_ax2 = fig.add_axes([0.606, 0.2, 0.015, 0.6])
        plt.colorbar(p2, cax=cbar_ax2, label="Rise time (s)")
        cbar_ax3 = fig.add_axes([0.893, 0.2, 0.015, 0.6])
        plt.colorbar(p3, cax=cbar_ax3, label="Rupture onset (s)")
        
        plt.subplots_adjust(wspace=0.4)
        
        name_out=home+project_name+'/plots/'+faults[kfault].split('/')[-1].split('.')[0]+'.'+faults[kfault].split('/')[-1].split('.')[1]+'.png'
        plt.savefig(name_out)
        plt.close()


def plot_LW_scaling(home,project_name,run_name):
    '''
    Plot a comparisson between used length/width and Blasser scaling laws
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10
    from matplotlib import pyplot as plt
    
    logs=glob(home+project_name+'/output/ruptures/'+run_name+'*.log')
    L=zeros(len(logs))
    W=zeros(len(logs))
    Mw=zeros(len(logs))
    for k in range(len(logs)):
        f=open(logs[k],'r')
        while True:
            line=f.readline()
            if 'Lmax' in line:
                L[k]=float(line.split(':')[-1].split()[0])
            if 'Wmax' in line:
                W[k]=float(line.split(':')[-1].split()[0])
            if 'Actual magnitude' in line:
                Mw[k]=float(line.split()[-1])
                break
    
    #Make plot
    plt.figure(figsize=(9,4))
    
    plt.subplot(121)
    Mw_synth=arange(7.8,9.3,0.1)
    plt.plot(Mw_synth,-2.37+0.57*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(L),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (L) [km]')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(8.0,3.17))
    
    plt.subplot(122)
    plt.plot(Mw_synth,-1.86+0.46*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(W),marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (W) [km]')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(8.0,2.8))
    
    plt.show()
    
    plt.subplots_adjust(bottom=0.17)
 
    
def analyze_sources(home,project_name,run_name,Mw_lims=[7.75,9.35]):
    '''
    Basic parameters of the sources
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10,genfromtxt,sqrt
    from matplotlib import pyplot as plt
    
    logs=glob(home+project_name+'/output/ruptures/'+run_name+'*.log')
    ruptures=glob(home+project_name+'/output/ruptures/'+run_name+'*.rupt')
    L=zeros(len(logs))
    W=zeros(len(logs))
    Mw_target=zeros(len(logs))
    Mw_actual=zeros(len(logs))
    slip_mean=zeros(len(logs))
    slip_max=zeros(len(logs))
    slip_stdev=zeros(len(logs))
    trise_mean=zeros(len(logs))
    trise_max=zeros(len(logs))
    trise_stdev=zeros(len(logs))
    for k in range(len(logs)):
        f=open(logs[k],'r')
        while True:
            line=f.readline()
            if 'Lmax' in line:
                L[k]=float(line.split(':')[-1].split()[0])
            if 'Wmax' in line:
                W[k]=float(line.split(':')[-1].split()[0])
            if 'Target magnitude' in line:
                Mw_target[k]=float(line.split()[-1])
            if 'Actual magnitude' in line:
                Mw_actual[k]=float(line.split()[-1])
                break
        #now go look in rupture files for slip aprameters
        source=genfromtxt(ruptures[k])
        slip=sqrt(source[:,8]**2+source[:,9]**2)
        trise=source[:,7]
        slip_mean[k]=slip.mean()
        slip_max[k]=slip.max()
        slip_stdev[k]=slip.std()
        trise_mean[k]=trise.mean()
        trise_max[k]=trise.max()
        trise_stdev[k]=trise.std()
    #Make plot
    plt.figure(figsize=(15,10.5))
    
    plt.subplot(331)
    Mw_synth=arange(7.8,9.3,0.1)
    plt.plot(Mw_synth,-2.37+0.57*Mw_synth,c='k',lw=2)
    plt.scatter(Mw_actual,log10(L),marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('log (L) [km]')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(7.9,3.0))
    
    plt.subplot(332)
    plt.plot(Mw_synth,-1.86+0.46*Mw_synth,c='k',lw=2)
    plt.scatter(Mw_actual,log10(W),marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('log (W) [km]')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(7.9,2.4))
    
    plt.subplot(333)
    plt.scatter(Mw_actual,Mw_target,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Target Mw')
    
    plt.subplot(334)
    plt.scatter(Mw_actual,slip_mean,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Mean slip (m)')
    
    plt.subplot(335)
    plt.scatter(Mw_actual,slip_max,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Peak slip (m)')
    
    plt.subplot(336)
    plt.scatter(Mw_actual,slip_stdev,marker='+')
    plt.xlim([7.75,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('Slip Std.Dev. (m)')
    
    plt.subplot(337)
    plt.scatter(Mw_actual,trise_mean,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Mean rise time (s)')
    
    plt.subplot(338)
    plt.scatter(Mw_actual,trise_max,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Peak rise time (s)')
    
    plt.subplot(339)
    plt.scatter(Mw_actual,trise_stdev,marker='+')
    plt.xlim(Mw_lims)
    plt.xlabel('Actual Mw')
    plt.ylabel('Rise time Std.Dev. (s)')
    
    plt.show()
    
    plt.subplots_adjust(bottom=0.1,left=0.06,right=0.97,top=0.97)
    
                        

def one_event_pgd_scaling(home,project_name,run_name,run_number,reference='centroid',dist_lims=[1,1000],plus_minus=[0.1,0.2,0.3]):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,log10
    from matplotlib import pyplot as plt
    from string import replace
    from mudpy.analysis import pgd_regression
    
    # Read summary file
    summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    #summary_file=summary_file=home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'

    lonlat=genfromtxt(summary_file,usecols=[1,2])
    pgd=genfromtxt(summary_file,usecols=[6])*100
    
    #Second set of pgd's
    #summary_file2=home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    summary_file2='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    pgd2=genfromtxt(summary_file2,usecols=[6])*100
    
    # Get hypocenter or centroid
    event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if reference=='hypocenter':
            if 'Hypocenter (lon,lat,z[km])' in line:
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False
        elif reference=='centroid':
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
        if 'Actual magnitude' in line:
            Mw=float(line.split(':')[-1].split(' ')[-1])

    #compute station to hypo distances
    d=zeros(len(lonlat))
    for k in range(len(lonlat)):
        d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
        d[k]=d[k]/1000
        
    # Get best fitting law for this event
    A,B,C=pgd_regression(home,project_name,run_name,run_number)
        
        
    fig = plt.figure()
    plt.title('Mw = '+str(Mw))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    # Plot reference line
    dref=logspace(0,3,50)
    pgdref=10**(-4.434+1.047*Mw-0.138*Mw*log10(dref))
    plt.plot(dref,pgdref,'k',lw=2,label='Melgar et al.')
    pgdref=10**(A+B*Mw+C*Mw*log10(dref))
    plt.plot(dref,pgdref,'r',lw=2,label='Event best fit')
    #Plot plus_minus lines
    for k in range(len(plus_minus)):
        Mw_plus=Mw+plus_minus[k]
        pgdref_plus=10**(-4.434+1.047*Mw_plus-0.138*Mw_plus*log10(dref))
        plt.plot(dref,pgdref_plus,'#A0A0A0',lw=1.0,label=None)
        Mw_minus=Mw-plus_minus[k]
        pgdref_minus=10**(-4.434+1.047*Mw_minus-0.138*Mw_minus*log10(dref))
        plt.plot(dref,pgdref_minus,'#A0A0A0',lw=1.5,label=None)
    #Actual observations
    ax.scatter(d,pgd,s=60, facecolors='none', edgecolors='m',marker='d',label='Simulation')
    #ax.scatter(d,pgd2,s=60, facecolors='none', edgecolors='k',marker='d',label='gil7.mod')
    plt.legend(loc=3)
    plt.xlabel('Station to event distance (km)')
    plt.ylabel('PGD (cm)')
    plt.xlim(dist_lims)
    plt.ylim([pgd.min(),pgd.max()+100])
    plt.subplots_adjust(bottom=0.15)
    plt.show()
   
    
     
def one_event_pgd_distance_misfit(home,project_name,run_name,run_number,dist_lims=[1,1000],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-2,2]):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,log10,log
    from matplotlib import pyplot as plt
    from string import replace
    from mudpy.analysis import pgd_regression
    
    # Read summary file
    summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    lonlat=genfromtxt(summary_file,usecols=[1,2])
    pgd=genfromtxt(summary_file,usecols=[6])*100
    
    # Get hypocenter or centroid
    event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'Centroid (lon,lat,z[km])' in line:                
            s=replace(line.split(':')[-1],'(','')
            s=replace(s,')','')
            hypo=array(s.split(',')).astype('float')
            loop_go=False       
        if 'Actual magnitude' in line:
            Mw=float(line.split(':')[-1].split(' ')[-1])

    #compute station to hypo distances
    d=zeros(len(lonlat))
    for k in range(len(lonlat)):
        d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
        d[k]=d[k]/1000
         
    # Model pgd           
    model_pgd=10**(A+B*Mw+C*Mw*log10(d))
        
        
    # Misfit
    misfit=log(pgd/model_pgd)
        
    fig = plt.figure()
    plt.title(run_name+'.'+run_number+' (Mw = '+str(Mw)+')')
    ax = plt.gca()
    ax.scatter(d,misfit,s=80, facecolors='none', edgecolors='k',marker='d')
    ax.set_xscale('log')
    # Plot reference line
    plt.xlabel('Station to event distance (km)')
    plt.ylabel('ln(Data/Model)')
    plt.xlim(dist_lims)

    plt.ylim(misfit_lims)
    plt.subplots_adjust(bottom=0.15)
    plt.show()   
  
        
              
                    
def pgd_distance_misfit(home,project_name,hist='True',rupture_list='ruptures.list',dist_lims=[10,1000],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-3,3],nbins=10):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,r_,log,genfromtxt,log10,ones,where
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from string import replace
    from matplotlib.ticker import MultipleLocator,LogLocator
    
    xmajorLocator = LogLocator(base=10.0, subs=[1])
    xminorLocator = MultipleLocator(1)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.2)
    
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    pgd_all=array([])
    d_all=array([])
    Mw_all=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
        
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all=r_[pgd_all,pgd]
        d_all=r_[d_all,d]
        Mw_all=r_[Mw_all,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]                         
    
    
    #Second set of pgd's
    pgd_all2=array([])
    d_all=array([])
    Mw_all=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        #summary_file=summary_file=home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        #summary_file='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        summary_file='/Users/dmelgar/FakeQuakes/Cascadia_nolognormal/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
        
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all2=r_[pgd_all2,pgd]
        d_all=r_[d_all,d]
        Mw_all=r_[Mw_all,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]    
    
    
    
    #get misfit
    misfit=log(pgd_all/pgd_predicted_all)
    misfit2=log(pgd_all2/pgd_predicted_all)
    
    if hist==False:                                                                           
        fig = plt.figure()
        ax = plt.gca()
        s=ax.scatter(d_all,misfit,c=Mw_all,cmap=cm.nipy_spectral_r,s=10,lw=0)
        cb=plt.colorbar(s)
        cb.set_label('Mw')
        ax.set_xscale('log')
        # Plot reference line
        plt.xlabel('Station to event distance (km)')
        plt.ylabel('ln(Data/Model)')
        plt.xlim(dist_lims)
    
        plt.ylim(misfit_lims)
        plt.subplots_adjust(bottom=0.15)
    else:
        bin_edges=logspace(log10(dist_lims[0]),log10(dist_lims[1]),nbins)
        bin_centers=zeros(len(bin_edges)-1)
        misfit_bin=zeros(len(bin_edges)-1)
        misfit_bin2=zeros(len(bin_edges)-1)
        misfit_stdev=zeros(len(bin_edges)-1)
        for k in range(len(bin_edges)-1):
            i=where((d_all>=bin_edges[k]) & (d_all<bin_edges[k+1]))[0]
            misfit_bin[k]=misfit[i].mean()
            misfit_bin2[k]=misfit2[i].mean()
            misfit_stdev[k]=misfit[i].std()
            bin_centers[k]=10**(log10(bin_edges[k+1])-(log10(bin_edges[k+1])-log10(bin_edges[k]))/2)
        
        fig = plt.figure(figsize=(9,5))
        ax = plt.gca()
        ax.plot(bin_centers,misfit_bin,c='k',lw=2,label='log-normal')
        #ax.plot(bin_centers,misfit_bin2,c='r',lw=2,label='normal')
        plt.legend(loc=2)
        ax.scatter(d_all,misfit,c='#909090',s=5,lw=0)
        ax.scatter(bin_centers,misfit_bin,c='k',lw=0.5,s=70)
        ax.plot(bin_centers,misfit_bin+1.64*misfit_stdev,c='k',lw=1)
        ax.plot(bin_centers,misfit_bin-1.64*misfit_stdev,c='k',lw=1)
        #ax.scatter(bin_centers,misfit_bin2,c='r',lw=0.5,s=70)
        ax.set_xscale('log')
        ax = plt.gca()
        ax.xaxis.set_major_locator(xmajorLocator)
        #ax.xaxis.set_minor_locator(xminorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)
        ax.tick_params(which='major',length=7,width=1)
        ax.tick_params(which='minor',length=4,width=1)
        # Plot reference line
        plt.xlabel('Station to event distance (km)')
        plt.ylabel('Ln(Simulation / GMPE)')
        plt.xlim(dist_lims)
        plt.ylim(misfit_lims)
        plt.subplots_adjust(bottom=0.15)
        plt.grid()
    plt.show()                                               
                                                        
 
 
def pgd_magnitude_misfit(home,project_name,hist=True,rupture_list='ruptures.list',Mw_lims=[7.8,9.5],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-3,3],nbins=10):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,linspace,r_,log,genfromtxt,log10,ones,where
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from string import replace
    from matplotlib.ticker import MultipleLocator
    
    xmajorLocator = MultipleLocator(0.5)
    xminorLocator = MultipleLocator(0.1)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.2)
    
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    pgd_all=array([])
    d_all=array([])
    Mw_all=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
                print Mw
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
        
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all=r_[pgd_all,pgd]
        d_all=r_[d_all,d]
        Mw_all=r_[Mw_all,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]                         
    
    
    #Second set of pgd's
    pgd_all2=array([])
    d_all=array([])
    Mw_all2=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        #summary_file=summary_file=home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        #summary_file='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        summary_file='/Users/dmelgar/FakeQuakes/Cascadia_nolognormal/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
        
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all2=r_[pgd_all2,pgd]
        d_all=r_[d_all,d]
        Mw_all2=r_[Mw_all2,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]    
    
    
    
    #get misfit
    misfit=log(pgd_all/pgd_predicted_all)
    misfit2=log(pgd_all2/pgd_predicted_all)
    
    if hist==False:                                                                           
        fig = plt.figure()
        ax = plt.gca()
        s=ax.scatter(d_all,misfit,c=Mw_all,cmap=cm.nipy_spectral_r,s=10,lw=0)
        cb=plt.colorbar(s)
        cb.set_label('Mw')
        ax.set_xscale('log')
        # Plot reference line
        plt.xlabel('Station to event distance (km)')
        plt.ylabel('ln(Data/Model)')
        plt.xlim(Mw_lims)
    
        plt.ylim(misfit_lims)
        plt.subplots_adjust(bottom=0.15)
    else:
        bin_edges=linspace(Mw_all.min(),Mw_all.max(),nbins)
        bin_centers=zeros(len(bin_edges)-1)
        misfit_bin=zeros(len(bin_edges)-1)
        misfit_stdev=zeros(len(bin_edges)-1)
        misfit_bin2=zeros(len(bin_edges)-1)
        for k in range(len(bin_edges)-1):
            i=where((Mw_all>=bin_edges[k]) & (Mw_all<bin_edges[k+1]))[0]
            misfit_bin[k]=misfit[i].mean()
            misfit_stdev[k]=misfit[i].std()
            misfit_bin2[k]=misfit2[i].mean()
            bin_centers[k]=bin_edges[k+1]-((bin_edges[k+1]-bin_edges[k])/2)
        
        fig = plt.figure(figsize=(9,5))
        ax = plt.gca()
        #Plot lines for legend
        ax.plot(bin_centers,misfit_bin,c='k',lw=2,label='log-normal')
        #ax.plot(bin_centers,misfit_bin2,c='r',lw=2,label='normal')
        plt.legend(loc=2)
        #Plot all events in background
        ax.scatter(Mw_all,misfit,c='#909090',s=5,lw=0)
        ax.scatter(bin_centers,misfit_bin,c='k',lw=0.5,s=70)
        ax.plot(bin_centers,misfit_bin+1.64*misfit_stdev,c='k',lw=1)
        ax.plot(bin_centers,misfit_bin-1.64*misfit_stdev,c='k',lw=1)
        #Plot again over the event ones
        ax.plot(bin_centers,misfit_bin,c='k',lw=2)
        #ax.plot(bin_centers,misfit_bin2,c='r',lw=2)
        ax.scatter(bin_centers,misfit_bin,c='k',lw=0.5,s=70)
        #ax.scatter(bin_centers,misfit_bin2,c='r',lw=0.5,s=70)
        # Plot reference line
        plt.xlabel('Event magnitude (Mw)')
        plt.ylabel('Ln(Simulation / GMPE)')
        plt.xlim(Mw_lims)
        plt.ylim(misfit_lims)
        plt.subplots_adjust(bottom=0.15)
        ax = plt.gca()
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_minor_locator(xminorLocator)
        ax.yaxis.set_major_locator(ymajorLocator)
        ax.yaxis.set_minor_locator(yminorLocator)
        ax.tick_params(which='major',length=7,width=1)
        ax.tick_params(which='minor',length=4,width=1)
        plt.grid()
    plt.show()        
   
         
def pgd_2D_misfit(home,project_name,rupture_list='ruptures.list',Mw_lims=[7.8,9.5],dist_lims=[10,1000],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-3,3],n_mag_bins=10,n_dist_bins=10):
    '''
    Plot misfit as a function of both distance and magnitude
    '''
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,linspace,r_,log,genfromtxt,log10,ones,where,arange,median,mean
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from string import replace
    from matplotlib.ticker import MultipleLocator,LogLocator,MaxNLocator
    
    xmajorLocator = LogLocator(base=10.0, subs=[1])
    xminorLocator = MultipleLocator(1)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.2)
    
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    pgd_all=array([])
    d_all=array([])
    Mw_all=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        # summary_file=summary_file=home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
           
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all=r_[pgd_all,pgd]
        d_all=r_[d_all,d]
        Mw_all=r_[Mw_all,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]               
    #Get misfits             
    misfit=log(pgd_all/pgd_predicted_all)                    
   
    #plotting
    bin_edges_x=linspace(Mw_all.min(),Mw_all.max(),n_mag_bins)
    bin_centers_x=zeros(len(bin_edges_x)-1)
    bin_edges_y=logspace(log10(dist_lims[0]),log10(dist_lims[1]),n_dist_bins)
    bin_centers_y=zeros(len(bin_edges_y)-1)
    misfit_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    stdev_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    frequency_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
        
    for kx in range(len(bin_centers_x)):
        for ky in range(len(bin_centers_y)):
            i=where((Mw_all>=bin_edges_x[kx]) & (Mw_all<bin_edges_x[kx+1]) & (d_all>=bin_edges_y[ky]) & (d_all<bin_edges_y[ky+1]))[0]
            #misfit_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit))
            misfit_bin[kx,ky]=mean(misfit[i])
            #misfit_bin[kx,ky]=median(misfit[i])
            stdev_bin[kx,ky]=misfit[i].std()
            frequency_bin[kx,ky]=len(misfit[i])
            bin_centers_x[kx]=bin_edges_x[kx+1]-((bin_edges_x[kx+1]-bin_edges_x[kx])/2)        
            bin_centers_y[ky]=bin_edges_y[ky+1]-((bin_edges_y[ky+1]-bin_edges_y[ky])/2) 
            
    fig=plt.figure(figsize=(15,5)) 
    
    plt.subplot(131)                           
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,misfit_bin.T,vmin=-1.5,vmax=1.5,cmap=cm.RdBu_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,misfit_bin.T,colors='#505050',levels=arange(-3,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    plt.ylabel('Distance (km)')
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('Ln(Simulation / GMPE)')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(132)
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,stdev_bin.T,vmin=0,vmax=1.5,cmap=cm.gist_heat_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,misfit_bin.T,colors='#505050',levels=arange(-3,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    #plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    plt.tick_params(axis='y',labelleft='off')
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('Misfit Std. Dev.')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=5)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(133)
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,frequency_bin.T,cmap=cm.magma_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,misfit_bin.T,colors='#505050',levels=arange(-3,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    #plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    plt.tick_params(axis='y',labelleft='off')
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('# of waveforms')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplots_adjust(left=0.08,top=0.88,right=0.97,bottom=0.1,wspace=0.07)
                                                                   

def pgd_2D_GOF(home,project_name,rupture_list='ruptures.list',Mw_lims=[7.8,9.5],dist_lims=[10,1000],A=-4.434,B=1.047,C=-0.138,GOF_lims=[0,2],n_mag_bins=10,n_dist_bins=10):
    '''
    Plot misfit as a function of both distance and magnitude
    '''
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,linspace,r_,log,genfromtxt,log10,ones,where,arange,median,mean
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from string import replace
    from matplotlib.ticker import MultipleLocator,LogLocator,MaxNLocator
    
    xmajorLocator = LogLocator(base=10.0, subs=[1])
    xminorLocator = MultipleLocator(1)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.2)
    
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    pgd_all=array([])
    d_all=array([])
    Mw_all=array([])
    pgd_predicted_all=array([])
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read summary file
        summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        # summary_file=summary_file=home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
    
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000   
           
        #Get predicted
        pgd_predicted=10**(A+B*Mw+C*Mw*log10(d))
        #Concatente to output variables
        pgd_all=r_[pgd_all,pgd]
        d_all=r_[d_all,d]
        Mw_all=r_[Mw_all,Mw*ones(len(d))]                   
        pgd_predicted_all=r_[pgd_predicted_all,pgd_predicted]               
    #Get misfits             
    misfit=log(pgd_all/pgd_predicted_all)                    
   
    #plotting
    bin_edges_x=linspace(Mw_all.min(),Mw_all.max(),n_mag_bins)
    bin_centers_x=zeros(len(bin_edges_x)-1)
    bin_edges_y=logspace(log10(dist_lims[0]),log10(dist_lims[1]),n_dist_bins)
    bin_centers_y=zeros(len(bin_edges_y)-1)
    misfit_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    stdev_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    frequency_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
        
    for kx in range(len(bin_centers_x)):
        for ky in range(len(bin_centers_y)):
            i=where((Mw_all>=bin_edges_x[kx]) & (Mw_all<bin_edges_x[kx+1]) & (d_all>=bin_edges_y[ky]) & (d_all<bin_edges_y[ky+1]))[0]
            misfit_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit))
            stdev_bin[kx,ky]=misfit[i].std()
            frequency_bin[kx,ky]=len(misfit[i])
            bin_centers_x[kx]=bin_edges_x[kx+1]-((bin_edges_x[kx+1]-bin_edges_x[kx])/2)        
            bin_centers_y[ky]=bin_edges_y[ky+1]-((bin_edges_y[ky+1]-bin_edges_y[ky])/2) 
            
    fig=plt.figure(figsize=(10,5)) 
    
    plt.subplot(121)                           
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,misfit_bin.T,vmin=GOF_lims[0],vmax=GOF_lims[1],cmap=cm.afmhot_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,misfit_bin.T,colors='#808080',levels=arange(-0,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    plt.ylabel('Distance (km)')
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('CGOF')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(122)
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,frequency_bin.T,cmap=cm.magma_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,misfit_bin.T,colors='#808080',levels=arange(0,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    #plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    plt.tick_params(axis='y',labelleft='off')
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('# of waveforms')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplots_adjust(left=0.08,top=0.88,right=0.97,bottom=0.1,wspace=0.07)           
                                 
                
                          
def record_section(home,project_name,GF_list,rupture,factor=10):
    '''
    Plot record section
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,ones,zeros
    from matplotlib import pyplot as plt
    from string import replace
    from obspy import read
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    xmajorLocator = MultipleLocator(50)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(200)
    yminorLocator = MultipleLocator(50)
    
    # Read summary file
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    event_log=home+project_name+'/output/ruptures/'+rupture+'.log'
    
    #Load velocity model for ray tracing
    velmod = velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
    
    # Get hypocenter
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'Hypocenter (lon,lat,z[km])' in line:
            s=replace(line.split(':')[-1],'(','')
            s=replace(s,')','')
            hypo=array(s.split(',')).astype('float')
            loop_go=False

    #compute station to hypo distances
    d=zeros(len(lonlat))
    for k in range(len(lonlat)):
        d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
        d[k]=d[k]/1000
        if lonlat[k,1]<hypo[1]: #station is south
            d[k]=-d[k]
 
    ptime=1e6*ones(len(sta))
    stime=1e6*ones(len(sta))
    for k in range(len(sta)):
        offset=d[k]
        e=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYE.sac')
        n=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYN.sac')
        z=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYZ.sac')            
        
        #Normalize and scale and apply distance offset
        n[0].data=factor*(n[0].data/max(abs(n[0].data)))+offset
        e[0].data=factor*(e[0].data/max(abs(e[0].data)))+offset
        z[0].data=factor*(z[0].data/max(abs(z[0].data)))+offset
        
        if k==0:
            N=n.copy()
            E=e.copy()
            Z=z.copy()
        else:
            N+=n
            E+=e
            Z+=z
        
        #Get ray arrivals
        deg=locations2degrees(hypo[1],hypo[0],lonlat[k,1],lonlat[k,0])
        
        # Ray trace
        try: #This is to deal with some weird bug in obspy.TauP
            arrivals = velmod.get_travel_times(source_depth_in_km=hypo[2],distance_in_degree=deg,phase_list=['P','Pn','S','Sn','p','s'])
        except:
            print 'No phases for '+rupture
            plt.close()
            return
        
        #Determine P and S arrivals
        for kphase in range(len(arrivals)):
            if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
                if arrivals[kphase].time<ptime[k]:
                    ptime[k]=arrivals[kphase].time
            if 'S' == arrivals[kphase].name or 's' == arrivals[kphase].name or 'Sn' == arrivals[kphase].name:
                if arrivals[kphase].time<stime[k]:
                    stime[k]=arrivals[kphase].time
    
    #plot on figure
    plt.figure(figsize=(18,7))
    
    plt.subplot(131)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(N[k].times(),N[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([n[0].times()[0],n[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    plt.ylabel('Hypocentral distance (km)')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.title('North')
    
    plt.subplot(132)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(E[k].times(),E[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([e[0].times()[0],e[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.tick_params(axis='y',labelleft='off')
    plt.title('East')
    
    plt.subplot(133)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(Z[k].times(),Z[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([z[0].times()[0],z[0].times()[-1]])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    plt.xlabel('Seconds since OT')
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.tick_params(axis='y',labelleft='off')
    plt.title('Vertical')
    
    plt.subplots_adjust(left=0.05,right=0.98,wspace=0.05)
    plt.show()
    
    
    
def statics_map(summary_file,scale,xl,yl,vertscale=400):
    '''
    Plot statics
    '''
    
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,zeros,r_,where
    
    sta=genfromtxt(summary_file)
    
    
    plt.figure(figsize=(5,12))
    
    #Plot coasts (also hard coded to cascadia for now)
    fcoast=open('/Users/dmelgar/Cascadia/coast/coast.txt','r')
    line=fcoast.readline()
    parsing_coast=True
    while parsing_coast==True:
        if '>' in line:
            complete_polygon=False
            count=0
            while complete_polygon==False:
                line=fcoast.readline()
                if '>' in line or line=='':
                    complete_polygon=True
                    plt.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.3,c='k')
                else:
                    new_coordinates=zeros((1,2))
                    new_coordinates[0,0]=float(line.split()[0])
                    new_coordinates[0,1]=float(line.split()[1])
                    if count==0:
                        coast_coordinates=new_coordinates.copy()
                        count+=1
                    else:
                        coast_coordinates=r_[coast_coordinates,new_coordinates]
                if line=='':
                    parsing_coast=False
    #Done plotting coast
    #Plot cnational boundaries(also hard coded to cascadia for now)
    fcoast=open('/Users/dmelgar/Cascadia/coast/boundaries.txt','r')
    line=fcoast.readline()
    parsing_coast=True
    while parsing_coast==True:
        if '>' in line:
            complete_polygon=False
            count=0
            while complete_polygon==False:
                line=fcoast.readline()
                if '>' in line or line=='':
                    complete_polygon=True
                    plt.plot(coast_coordinates[:,0],coast_coordinates[:,1],lw=0.5,c='k')
                else:
                    new_coordinates=zeros((1,2))
                    new_coordinates[0,0]=float(line.split()[0])
                    new_coordinates[0,1]=float(line.split()[1])
                    if count==0:
                        coast_coordinates=new_coordinates.copy()
                        count+=1
                    else:
                        coast_coordinates=r_[coast_coordinates,new_coordinates]
                if line=='':
                    parsing_coast=False
    #Done plotting coast
    
    iup=where(sta[:,5]>0)
    idown=where(sta[:,5]<0)
    plt.scatter(sta[iup,1],sta[iup,2],c='r',s=sta[iup,5]*vertscale,lw=0.1)
    plt.scatter(sta[idown,1],sta[idown,2],c='#4169E1',s=-sta[idown,5]*vertscale,lw=0.1)
    plt.quiver(sta[:,1],sta[:,2],sta[:,4],sta[:,3],scale=scale,width=0.007)
    
    plt.xlim(xl)
    plt.ylim(yl)
    
    #Reference
    ref=1
    plt.quiver(-123.8,50.7,ref,0,scale=scale,width=0.007)
    plt.annotate('1m',xy=(-121.5,50.6))
    
    ref=1
    plt.scatter(-122.6,50.1,c='r',s=ref*vertscale,lw=0.1)
    plt.quiver(-122.6,49.98,0,1,10,width=0.007)
    plt.scatter(-123.6,50.1,c='#4169E1',s=ref*vertscale,lw=0.1)
    plt.quiver(-123.6,50.25,0,-1,10,width=0.007)
    plt.annotate('1m',xy=(-121.5,50.0))
    
    plt.show()
    


def plot_PGD_all(home,project_name,rupture_list,dist_range=[1,1000],Mw_range=[6,9.5]):
    '''
    Plot comaprisson between scaling law and observed PGD's
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from numpy import genfromtxt,zeros,array,r_,ones,log10,arange,concatenate
    from string import replace
    from matplotlib.collections import LineCollection
    
    #get list of ruptures
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    #Init
    Mw_plot=array([])
    dist_plot=array([])
    pgd_plot=array([])
    for krupt in range(len(ruptures)):
        print 'Reading data from rupture '+ruptures[krupt]
        # Read summary file
        summary_file=home+project_name+'/output/waveforms/'+ruptures[krupt].split('.')[0]+'.'+ruptures[krupt].split('.')[1]+'/_summary.'+ruptures[krupt].split('.')[0]+'.'+ruptures[krupt].split('.')[1]+'.txt'
        lonlat=genfromtxt(summary_file,usecols=[1,2])
        pgd=genfromtxt(summary_file,usecols=[6])*100
        # Get hypocenter or centroid
        event_log=home+project_name+'/output/ruptures/'+ruptures[krupt].split('.')[0]+'.'+ruptures[krupt].split('.')[1]+'.log'
        f=open(event_log,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Centroid (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False       
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
        #compute station to hypo distances
        d=zeros(len(lonlat))
        for k in range(len(lonlat)):
            d[k],az,baz=gps2dist_azimuth(lonlat[k,1],lonlat[k,0],hypo[1],hypo[0])
            d[k]=d[k]/1000
        #Append to array of magnitudes, pgd's and distances
        Mw_plot=r_[Mw_plot,Mw*ones(len(d))]
        pgd_plot=r_[pgd_plot,pgd]
        dist_plot=r_[dist_plot,d]
    
    
    plt.figure()
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    s=ax.scatter(dist_plot,pgd_plot,c=Mw_plot,marker='o',s=5,lw=0,cmap=cm.nipy_spectral_r)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(18)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(18)
    cb=plt.colorbar(s)
    cb.set_label(r'$M_w$',fontsize=18)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(18)
        
    #Plot reference lines
#    Mw_ref=arange(Mw_range[0],Mw_range[-1]+0.1,0.1)
#    dist_ref=arange(dist_range[0],dist_range[-1])
#    for k in range(len(Mw_ref)):
#        pgd_ref=10**(-4.434+1.047*Mw_ref[k]-0.138*Mw_ref[k]*log10(dist_ref))
#        Mw_ref_plot=ones(len(dist_ref))*Mw_ref[k]
#        ##
#        points = array([dist_ref, pgd_ref]).T.reshape(-1, 1, 2)
#        segments = concatenate([points[:-1], points[1:]], axis=1)
#
#        lc = LineCollection(segments, cmap=plt.get_cmap('nipy_spectral_r'),)
#        lc.set_array(Mw_ref_plot)
#        lc.set_linewidth(3)
        
    
    plt.xlim(dist_range)
    plt.ylim([pgd_plot.min(),pgd_plot.max()])
    
    
def plot_pgd_decay(home,project_name,rupture_list='ruptures.list'):
    '''
    Plot observed decay vs predicted decay
    '''
    
    from mudpy.analysis import pgd_slope_difference
    from matplotlib import pyplot as plt
    
    Mw,theory,obs=pgd_slope_difference(home,project_name,rupture_list)
    
    plt.plot
    plt.scatter(Mw,obs,c='r',s=20,marker='s')
    plt.scatter(Mw,theory,c='k',s=20,marker='s')
    plt.legend(['Best fit','Melgar et al.'])
    plt.xlabel('Mw')
    plt.ylabel('log(PGD) distance decay')
    plt.show()
    

def amplitude_all(home,project_name,nbins_ampl=20,nbins_freq=20,ampl_lims=[1e-10,1e0],freq_lims=[1e5,0.5]):
    '''
    Plot all PSD's
    '''
    
    from glob import glob
    from numpy import load,logspace,log10,zeros,histogram2d
    from matplotlib import cm
    from matplotlib import pyplot as plt
    from obspy.imaging.cm import pqlx

    
    paths=glob(home+project_name+'/analysis/frequency/*')
    #Define bins
    bin_edges_x=logspace(log10(freq_lims[0]),log10(freq_lims[1]),nbins_freq)
    bin_edges_y=logspace(log10(ampl_lims[0]),log10(ampl_lims[1]),nbins_ampl)
    bin_centers_x=zeros(len(bin_edges_x)-1)
    bin_centers_y=zeros(len(bin_edges_y)-1)
    hit_count=zeros((len(bin_centers_x),len(bin_centers_y)))   
    for krupt in range(len(paths)):
        print paths[krupt]
        PSDs=glob(paths[krupt]+'/*.npz')
        for kpsd in range(len(PSDs)):
            psd=load(PSDs[kpsd])
            nspec=psd['npsd']**0.5
            espec=psd['epsd']**0.5
            uspec=psd['upsd']**0.5
            #Normalize
            nspec=nspec/nspec.max()
            espec=espec/espec.max()
            uspec=uspec/uspec.max()
            #Frequencies
            fn=psd['fn']
            fe=psd['fe']
            fu=psd['fu']
            #Get hit count
            H,xb,yb=histogram2d(fn,nspec, bins=(bin_edges_x,bin_edges_y))
            hit_count=hit_count+H
            
    plt.figure()
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,hit_count/hit_count.sum(),cmap=pqlx)
    ax = plt.gca()    
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude (m/Hz)^0.5')
    cb = plt.colorbar(PM) 
    cb.set_label('Probability')
            
            
        