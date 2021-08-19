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
        print(faults[kfault])
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
        
        name_out=home+project_name+'/plots/'+faults[kfault].split('/')[-1].split('.')[0]+'.'+faults[kfault].split('/')[-1].split('.')[1]+'.jpg'
        plt.savefig(name_out)
        plt.close()


def plot_LW_scaling(home,project_name,run_name):
    '''
    Plot a comparisson between used length/width and Blasser scaling laws
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10,sort
    from matplotlib import pyplot as plt
    
    logs=sort(glob(home+project_name+'/output/ruptures/'+run_name+'*.log'))
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
    Mw_synth=arange(7.0,9.3,0.1)
    plt.plot(Mw_synth,-2.37+0.57*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(L),marker='+')
    plt.xlim([7.3,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (L) [km]')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(8.0,3.17))
    
    plt.subplot(122)
    plt.plot(Mw_synth,-1.86+0.46*Mw_synth,c='k',lw=2)
    plt.scatter(Mw,log10(W),marker='+')
    plt.xlim([7.3,9.35])
    plt.xlabel('Actual Mw')
    plt.ylabel('log (W) [km]')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(8.0,2.8))
    
    plt.show()
    
    plt.subplots_adjust(bottom=0.17)
    
    return L,W,Mw
 
    
def analyze_sources(rupture_folder,Mw_lims=[7.75,9.35],return_values=False):
    '''
    Basic parameters of the sources
    '''
    
    from glob import glob
    from numpy import zeros,arange,log10,genfromtxt,sqrt,where,sort
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator
    
    plt.rcParams.update({'font.size': 12})
    logs=sort(glob(rupture_folder+'*.log'))
    ruptures=sort(glob(rupture_folder+'*.rupt'))

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
        i=where(slip>0)[0]
        slip_mean[k]=slip[i].mean()
        slip_max[k]=slip[i].max()
        slip_stdev[k]=slip.std()
        trise_mean[k]=trise[i].mean()
        trise_max[k]=trise[i].max()
        trise_stdev[k]=trise[i].std()
    #Make plot
    plt.figure(figsize=(15,10.5))
    
    plt.subplot(331)
    Mw_synth=arange(7.8,9.3,0.1)
    # plt.plot(Mw_synth,10**(-2.37+0.57*Mw_synth),c='k',lw=2)
    plt.scatter(Mw_actual,L,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([70,1100])
    plt.ylabel('Fault length (km)')
    plt.annotate(r'$\log (L)=-2.37+0.57M_w$',xy=(7.95,85),fontsize=12)
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(a)',xy=(7.9,765),fontsize=12)
    ax=plt.gca()
    ax.set_yscale('log')
    xmajorLocator = MultipleLocator(0.2)
    xminorLocator = MultipleLocator(0.05)
    ymajorLocator = MultipleLocator(0.2)
    yminorLocator = MultipleLocator(0.05)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)    
    
    
    plt.subplot(332)
    # plt.plot(Mw_synth,10**(-1.86+0.46*Mw_synth),c='k',lw=2)
    plt.scatter(Mw_actual,W,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([9,300])
    plt.ylabel('Fault width (km)')
    plt.annotate(r'$\log (W)=-1.86+0.46M_w$',xy=(7.9,12.2))
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(b)',xy=(7.9,190),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(0.2)
    yminorLocator = MultipleLocator(0.05)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)  
    
    plt.subplot(333)
    plt.scatter(Mw_actual,Mw_target,marker='+')
    plt.xlim(Mw_lims)
    plt.ylabel('Target Mw')
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(c)',xy=(7.9,8.8),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(0.2)
    yminorLocator = MultipleLocator(0.05)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
  
    plt.subplot(334)
    plt.scatter(Mw_actual,slip_mean,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([0.04,115])
    plt.ylabel(r'Mean slip (m)')
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(d)',xy=(7.9,45),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(5)
    yminorLocator = MultipleLocator(1)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
    #Allen Hayes scaling lines
    # plt.plot([7.0,9.6],[0.4,21.1],c='k',lw=2)
    # plt.plot([7.0,9.6],[0.4,7.87],c='r',lw=2)
    # plt.plot([7.0,9.6],[0.55,10.71],c='orange',lw=2)
    # plt.plot([7.0,9.6],[1.20,23.33],c='g',lw=2)
    # plt.plot([7.0,9.27],[1.08,39.16],c='violet',lw=2)
        
    plt.subplot(335)
    plt.scatter(Mw_actual,slip_max,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([0.7,110])
    plt.ylabel(r'Max. slip (m)')
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(e)',xy=(7.9,51),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(10)
    yminorLocator = MultipleLocator(2)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
    #Allen Hayes scaling lines
    # plt.plot([7.0,9.6],[1.07,75],c='k',lw=2)
    # plt.plot([7.0,9.6],[3.21,225],c='k',lw=2)
    
    plt.subplot(336)
    plt.scatter(Mw_actual,slip_stdev,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([0.2,20])
    plt.ylabel(r'Slip std. dev. (m)')
    plt.tick_params(axis='x',labelbottom='off')
    plt.annotate('(f)',xy=(7.9,12),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(2)
    yminorLocator = MultipleLocator(0.5)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
        
    plt.subplot(337)
    #plt.plot(Mw_synth,(4.308e-7)*(10**(1.5*Mw_synth+9.1))**(1./3),c='k',lw=2)
    plt.scatter(Mw_actual,trise_mean,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([0.9,30])
    plt.xlabel('Actual Mw')
    plt.ylabel(r'Mean rise time (s)')
    plt.annotate('(g)',xy=(7.9,17.3),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(5)
    yminorLocator = MultipleLocator(1)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
    
    plt.subplot(338)
    plt.scatter(Mw_actual,trise_max,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([3,115])
    plt.xlabel('Actual Mw')
    plt.ylabel(r'Max rise time (s)')
    plt.annotate('(h)',xy=(7.9,70),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(10)
    yminorLocator = MultipleLocator(2)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
    
    plt.subplot(339)
    plt.scatter(Mw_actual,trise_stdev,marker='+')
    plt.xlim(Mw_lims)
    plt.ylim([0.4,15])
    plt.xlabel('Actual Mw')
    plt.ylabel(r'Rise time std. dev (s)')
    plt.annotate('(i)',xy=(7.9,10.),fontsize=12)
    ax=plt.gca()
    ymajorLocator = MultipleLocator(2)
    yminorLocator = MultipleLocator(0.5)
    ax.set_yscale('log')
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    #ax.yaxis.set_major_locator(ymajorLocator)
    #ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1) 
    
    plt.subplots_adjust(bottom=0.1,left=0.06,right=0.99,top=0.97,hspace=0.03)    
    plt.show()
    
    if return_values==True:
        return Mw_actual,L,slip_mean
    
    
                        

def one_event_pgd_scaling(home,project_name,run_name,run_number,reference='centroid',dist_lims=[1,1000],
                coeffs=[-4.434,1.047,-0.138],plus_minus=[0.1,0.2,0.3],title_mag=None):
    '''
    Make PGD scaling plot
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,log10
    from matplotlib import pyplot as plt
    from mudpy.analysis import pgd_regression
    
    # Read summary file
    summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    #summary_file=summary_file=home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'

    lonlat=genfromtxt(summary_file,usecols=[1,2])
    pgd=genfromtxt(summary_file,usecols=[6])*100
    
    #Second set of pgd's
    #summary_file2=home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    #summary_file2='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    #pgd2=genfromtxt(summary_file2,usecols=[6])*100
    
    # Get hypocenter or centroid
    event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if reference=='hypocenter':
            if 'Hypocenter (lon,lat,z[km])' in line:
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
                hypo=array(s.split(',')).astype('float')
                loop_go=False
        elif reference=='centroid':
            if 'Centroid (lon,lat,z[km])' in line:                
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
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
        
        
    fig = plt.figure(figsize=(4.5,4.5))
    plt.title('Mw = '+str(Mw))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    # Plot reference line
    dref=logspace(0,3.3,50)
    A=coeffs[0]
    B=coeffs[1]
    C=coeffs[2]
    pgdref=10**(A+B*Mw+C*Mw*log10(dref))
    plt.plot(dref,pgdref,'k',lw=2,label='PGD GMPE')
    pgdref=10**(A+B*Mw+C*Mw*log10(dref))
    #plt.plot(dref,pgdref,'r',lw=2,label='Event best fit')
    #Plot plus_minus lines
    for k in range(len(plus_minus)):
        Mw_plus=Mw+plus_minus[k]
        pgdref_plus=10**(A+B*Mw_plus+C*Mw_plus*log10(dref))
        plt.plot(dref,pgdref_plus,'#A0A0A0',lw=1.0,label=None)
        Mw_minus=Mw-plus_minus[k]
        pgdref_minus=10**(A+B*Mw_minus+C*Mw_minus*log10(dref))
        plt.plot(dref,pgdref_minus,'#A0A0A0',lw=1.5,label=None)
    #Actual observations
    ax.scatter(d,pgd,s=60, facecolors='none', edgecolors='#1E90FF',marker='d',label='Simulation')
    #ax.scatter(d,pgd2,s=60, facecolors='none', edgecolors='k',marker='d',label='gil7.mod')
    plt.legend(loc=3,frameon=False)
    plt.xlabel('Station to event distance (km)')
    plt.ylabel('PGD (cm)')
    plt.xlim(dist_lims)
    plt.ylim([pgd.min(),pgd.max()+100])
    plt.subplots_adjust(bottom=0.15,left=0.17,right=0.99)
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
                print(Mw)
    
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
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
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
            misfit_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit[i]))
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
                                 

def pgd_2D_misfit_and_GOF(home,project_name,rupture_list='ruptures.list',Mw_lims=[7.8,9.5],dist_lims=[10,1000],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-3,3],GOF_lims=[0,2],n_mag_bins=10,n_dist_bins=10):
    '''
    Plot misfit as a function of both distance and magnitude
    '''
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,zeros,logspace,linspace,r_,log,genfromtxt,log10,ones,where,arange,median,mean
    from matplotlib import pyplot as plt
    from matplotlib import cm
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
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
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

    #remove extrneous magnitudes
    i=where((Mw_all>=Mw_lims[0]) & (Mw_all<=Mw_lims[1]))[0]
    Mw_all=Mw_all[i]
    misfit=misfit[i]
    d_all=d_all[i]
    #plotting
    bin_edges_x=linspace(Mw_all.min(),Mw_all.max(),n_mag_bins)
    bin_centers_x=zeros(len(bin_edges_x)-1)
    bin_edges_y=logspace(log10(dist_lims[0]),log10(dist_lims[1]),n_dist_bins)
    bin_centers_y=zeros(len(bin_edges_y)-1)
    misfit_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    GOF_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
    frequency_bin=zeros((len(bin_edges_x)-1,len(bin_edges_y)-1))
        
    for kx in range(len(bin_centers_x)):
        for ky in range(len(bin_centers_y)):
            i=where((Mw_all>=bin_edges_x[kx]) & (Mw_all<bin_edges_x[kx+1]) & (d_all>=bin_edges_y[ky]) & (d_all<bin_edges_y[ky+1]))[0]
            GOF_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit[i]))
            misfit_bin[kx,ky]=misfit[i].mean()
            frequency_bin[kx,ky]=len(misfit[i])
            bin_centers_x[kx]=bin_edges_x[kx+1]-((bin_edges_x[kx+1]-bin_edges_x[kx])/2)        
            bin_centers_y[ky]=bin_edges_y[ky+1]-((bin_edges_y[ky+1]-bin_edges_y[ky])/2) 
            
    fig=plt.figure(figsize=(15,5.5)) 
    
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
    levs=ax.get_xticks()
    ax.set_xticklabels(levs,rotation=-55)
    #plt.annotate('(a)',xy=(7.9,680),fontsize=16)
    bbox_props = dict(boxstyle="round", fc="w")
    ax.text(7.9, 350, "(a)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('Ln(Simulation / GMPE)')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(132)                           
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,GOF_bin.T,vmin=GOF_lims[0],vmax=GOF_lims[1],cmap=cm.afmhot_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,GOF_bin.T,colors='#808080',levels=arange(-0,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    levs=ax.get_xticks()
    ax.set_xticklabels(levs,rotation=-55)
    plt.xlabel('Magnitude')
    plt.tick_params(axis='y',labelleft='off')
    bbox_props = dict(boxstyle="round", fc="w")
    ax.text(7.9, 350, "(b)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('CGOF')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(133)
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,frequency_bin.T,cmap=cm.magma_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,GOF_bin.T,colors='#808080',levels=arange(0,3.1,0.5))
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    #plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    levs=ax.get_xticks()
    ax.set_xticklabels(levs,rotation=-55)
    plt.xlabel('Magnitude')
    plt.tick_params(axis='y',labelleft='off')
    bbox_props = dict(boxstyle="round", fc="w")
    ax.text(7.9, 350, "(c)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('# of waveforms')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplots_adjust(left=0.07,top=0.87,right=0.98,bottom=0.1,wspace=0.06)           
                           
 

def plot_2_waveforms_2_spectra(home,project_name,run_name,run_number,sta1,sta2,component='N',fmin=0.03,ref1offset=[-1.65,-1.65],ref2offset=[-3,-3]):
    '''
    Comapre spectral decay of several STFs
    
    viewFQ.plot_2_waveforsm_2_spectra(home,project_name,run_name,'001245','ONAB','BEND',component='N',fmin=0.1,ref1offset=[-1.6,-2.2],ref2offset=[-3,-3.7])

 ref1offset=[-4.8,-4.0],ref2offset=[-8,-6.3]

       
    '''                 
    
    from numpy import genfromtxt,r_,where,mean,log10,linspace
    from obspy import read
    from matplotlib import pyplot as plt
    from obspy.core import UTCDateTime
    import nitime.algorithms as tsa
    from matplotlib.ticker import MultipleLocator
    from scipy.signal import periodogram
    
    if component=='E':
        chan='LYE'
    elif component=='N':
        chan='LYN'
    elif component=='Z':
        chan='LYZ'
    else:
        print('Unknown channel')
        
    
    st1_triangle=read(home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/'+sta1+'.'+chan+'.sac')
    st1_cosine=read(home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/'+sta1+'.'+chan+'.sac')
    st1_dreger=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta1+'.'+chan+'.sac')
    #Nitime
    #f, psd1_triangle, nu = tsa.multi_taper_psd(st1_triangle[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #f, psd1_cosine, nu = tsa.multi_taper_psd(st1_cosine[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #f, psd1_dreger, nu = tsa.multi_taper_psd(st1_dreger[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #Periodogram
    f,psd1_triangle=periodogram(st1_triangle[0].data, fs=1./st1_triangle[0].stats.delta, nfft=None, return_onesided=True,window='hamming',scaling='spectrum')
    f,psd1_dreger=periodogram(st1_dreger[0].data, fs=1./st1_dreger[0].stats.delta, nfft=None, return_onesided=True,window='hamming',scaling='spectrum')
    #psd1_triangle=psd1_triangle**0.5
    #psd1_dreger=psd1_dreger**0.5
    
    st2_triangle=read(home+project_name+'/output/triangle/'+run_name+'.'+run_number+'/'+sta2+'.'+chan+'.sac')
    st2_cosine=read(home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/'+sta2+'.'+chan+'.sac')
    st2_dreger=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta2+'.'+chan+'.sac')
    #
    #f, psd2_triangle, nu = tsa.multi_taper_psd(st2_triangle[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #f, psd2_cosine, nu = tsa.multi_taper_psd(st2_cosine[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #f, psd2_dreger, nu = tsa.multi_taper_psd(st2_dreger[0].data,Fs=1./st1_triangle[0].stats.delta,adaptive=False,jackknife=False,low_bias=True)
    #
    f,psd2_triangle=periodogram(st2_triangle[0].data, fs=1./st2_triangle[0].stats.delta, nfft=None, return_onesided=True,window='hamming',scaling='spectrum')
    f,psd2_dreger=periodogram(st2_dreger[0].data, fs=1./st2_dreger[0].stats.delta, nfft=None, return_onesided=True,window='hamming',scaling='spectrum')
    #psd2_triangle=psd2_triangle**0.5
    #psd2_dreger=psd2_dreger**0.5    
    
    plt.figure(figsize=(10,6))
    
    ax=plt.axes([0.1,0.55,0.5,0.34])
    ax.xaxis.set_ticks_position('both')
    #ax.plot(st1_cosine[0].times(),st1_cosine[0].data,c='#303030',lw=1)
    ax.plot(st1_dreger[0].times(),st1_dreger[0].data,c='#606060',lw=1.5,label='Dreger STF')
    ax.plot(st1_triangle[0].times(),st1_triangle[0].data,c='#DC143C',lw=1.5,label='Triangle STF')
    plt.legend(frameon=False,ncol=2,labelspacing=0.1,bbox_to_anchor=(1.05, 1.25))
    ax.set_xlim([0,500])
    #ax.set_ylim([-0.3,0.3])
    #ax.yaxis.set_ticklabels(['','-0.1','','0','','0.1',''])
    ax.set_ylabel('North (m)')
    plt.tick_params(axis='x',labelbottom='off')
    xmajorLocator = MultipleLocator(100)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.1)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.annotate('Station '+sta1+'\nd=202km',xy=(300,1.7))
    
    ax=plt.axes([0.62,0.55,0.3,0.34])
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.tick_top()
    ax.xaxis.set_ticks_position('both')
    #plot amplitude spectrum
    ax.loglog(f,psd1_dreger,lw=1,c='#303030')
    ax.loglog(f,psd1_triangle,lw=1,c='#DC143C')
    #plot spectral ratio
    #ax.semilogx(f,(psd1_dreger**0.5)/(psd1_triangle**0.5),lw=1,c='#303030')
    #Plot reference lines
    fref=linspace(fmin,0.5)
    ref1=10**(-2*log10(fref)+ref1offset[0])
    ref2=10**(-3*log10(fref)+ref2offset[0])
    ax.loglog(fref,ref1,'--',c='#1E90FF',lw=3)
    ax.loglog(fref,ref2,'--',c='#1E90FF',lw=3)
    ax.set_xlim([0.002,0.5])
    ax.set_ylim([1e-8,5e-2])
    plt.tick_params(axis='x',labelbottom='off',labeltop='off')
    ax.set_ylabel('Amp. spec. (m/Hz)')
    ax.yaxis.set_label_position("right")
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    #Annotate reference lines
    ax.annotate(r'$f^{-2}$', xy=(0.025, 0.003))
    ax.annotate(r'$f^{-3}$', xy=(0.025, 0.00006))
    
    
    
    ax=plt.axes([0.1,0.15,0.5,0.34])
    #ax.plot(st2_cosine[0].times(),st2_cosine[0].data,c='#303030',lw=1)
    ax.plot(st2_dreger[0].times(),st2_dreger[0].data,c='#606060',lw=1.5)
    ax.plot(st2_triangle[0].times(),st2_triangle[0].data,c='#DC143C',lw=1.5)
    ax.set_xlim([0,500])
    #ax.set_ylim([-0.12,0.12])
    #ax.yaxis.set_ticklabels(['','-0.1','','0','','0.1',''])
    ax.set_ylabel('North (m)')
    ax.set_xlabel('Seconds after OT')
    xmajorLocator = MultipleLocator(100)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(0.2)
    yminorLocator = MultipleLocator(0.05)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.annotate('Station '+sta2+'\nd=383km',xy=(300,0.33))
    
    ax=plt.axes([0.62,0.15,0.3,0.34])
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    #plot spectra
    ax.loglog(f,psd2_dreger**0.5,lw=1,c='#303030')
    ax.loglog(f,psd2_triangle**0.5,lw=1,c='#DC143C')
    #plot spcertral ratios
    #ax.semilogx(f,(psd2_dreger**0.5)/(psd2_triangle**0.5),lw=1,c='#303030')
    #Plot reference lines
    fref=linspace(fmin,0.5)
    ref1=10**(-2*log10(fref)+ref1offset[1])
    ref2=10**(-3*log10(fref)+ref2offset[1])
    ax.loglog(fref,ref1,'--',c='#1E90FF',lw=3)
    ax.loglog(fref,ref2,'--',c='#1E90FF',lw=3)
    ax.set_xlim([0.002,0.5])
    ax.set_ylim([1e-6,5e-2])
    ax.set_ylabel('Amp. spec. (m/Hz)')
    ax.yaxis.set_label_position("right")
    ax.set_xlabel('Frequency (Hz)')
    #ax.yaxis.set_minor_locator(yminorLocator)      
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)        
                                                                
                          
def record_section(home,project_name,GF_list,rupture,factor=10,max_time=300):
    '''
    Plot record section
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,ones,zeros
    from matplotlib import pyplot as plt
    from obspy import read
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    
    xmajorLocator = MultipleLocator(100)
    xminorLocator = MultipleLocator(20)
    ymajorLocator = MultipleLocator(200)
    yminorLocator = MultipleLocator(50)
    
    # Read summary file
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    event_log=home+project_name+'/output/ruptures/'+rupture+'.log'
    
    #Load velocity model for ray tracing
#    velmod = TauPyModel(model="/Users/dmelgarm/FakeQuakes/Cascadia/structure/cascadia")
    velmod = TauPyModel(model="PREM")
    
    # Get hypocenter
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'Hypocenter (lon,lat,z[km])' in line:
            s=line.split(':')[-1]
            s=s.replace('(','')
            s=s.replace(')','')
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
        
        ## Ray trace
        #try: #This is to deal with some weird bug in obspy.TauP
        #    arrivals = velmod.get_travel_times(source_depth_in_km=hypo[2],distance_in_degree=deg,phase_list=['P','Pn','S','Sn','p','s'])
        #except:
        #    print('No phases for '+rupture
        #    plt.close()
        #    return
        
        ##Determine P and S arrivals
        #for kphase in range(len(arrivals)):
        #    if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
        #        if arrivals[kphase].time<ptime[k]:
        #            ptime[k]=arrivals[kphase].time
        #    if 'S' == arrivals[kphase].name or 's' == arrivals[kphase].name or 'Sn' == arrivals[kphase].name:
        #        if arrivals[kphase].time<stime[k]:
        #            stime[k]=arrivals[kphase].time
    
    #plot on figure
    plt.figure(figsize=(18,10))
    
    plt.subplot(131)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(N[k].times(),N[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([n[0].times()[0],max_time])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
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
    plt.xlim([e[0].times()[0],max_time])
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
    plt.yticks([])
    
    plt.subplot(133)
    plt.scatter(ptime,d,c='b',marker='|',s=80,lw=1.5)     
    plt.scatter(stime,d,c='r',marker='|',s=80,lw=1.5)    
    for k in range(len(sta)):
        plt.plot(Z[k].times(),Z[k].data,'k',lw=0.5)         
    #After axis adjsutments
    plt.xlim([z[0].times()[0],max_time])
    d_adjust=(d.max()-d.min())*0.05
    plt.ylim([d.min()-d_adjust,d.max()+d_adjust])
    ax = plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.tick_params(axis='y',labelleft='off')
    plt.title('Vertical')
    plt.yticks([])
    
    plt.subplots_adjust(left=0.15,right=0.98,bottom=0.07,top=0.96,wspace=0.05)
    plt.show()
    
    
    
def make_rgb_image(home,project_name,GF_list,rupture,factor=10,tlims=[0,300],saturation=None,lower_bound=None):
    '''
    Plot RGB image
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,ones,zeros,argsort,where,log10
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
    lon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1])
    lat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[2])
    event_log=home+project_name+'/output/ruptures/'+rupture+'.log'
    
    #sort stations
    i=argsort(lat)
    lat=lat[i]
    lon=lon[i]
    sta=sta[i]
    
    #Load velocity model for ray tracing
    velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
    
    # Get hypocenter
    f=open(event_log,'r')
    loop_go=True
    while loop_go:
        line=f.readline()
        if 'Actual magnitude' in line:
            Mw=float(line.split()[-1])
        if 'Hypocenter (lon,lat,z[km])' in line:
            s=replace(line.split(':')[-1],'(','')
            s=replace(s,')','')
            hypo=array(s.split(',')).astype('float')
            loop_go=False

    #compute station to hypo distances
    d=zeros(len(lon))
    for k in range(len(lon)):
        d[k],az,baz=gps2dist_azimuth(lat[k],lon[k],hypo[1],hypo[0])
        d[k]=d[k]/1000
        if lat[k]<hypo[1]: #station is south
            d[k]=-d[k]
          
    # sort by distance to hypo  
    i=argsort(d)
    lat=lat[i]
    lon=lon[i]
    sta=sta[i]    
 
    ptime=1e6*ones(len(sta))
    stime=1e6*ones(len(sta))
    pgd_east=0
    pgd_north=0
    pgd_up=0
    for k in range(len(sta)):

        e=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYE.sac')
        n=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYN.sac')
        z=read(home+project_name+'/output/waveforms/'+rupture+'/'+sta[k]+'.LYZ.sac')            
        
        # Beware this will not work if plotting based on sign
        if saturation!=None:
            i=where(abs(n[0].data)>saturation)[0]
            n[0].data[i]=saturation
            i=where(abs(e[0].data)>saturation)[0]
            e[0].data[i]=saturation
            i=where(abs(z[0].data)>saturation)[0]
            z[0].data[i]=saturation
            
        if lower_bound!=None:
            i=where(abs(n[0].data)<lower_bound)[0]
            n[0].data[i]=lower_bound
            i=where(abs(e[0].data)<lower_bound)[0]
            e[0].data[i]=lower_bound
            i=where(abs(z[0].data)<lower_bound)[0]
            z[0].data[i]=lower_bound

        if max(abs(n[0].data))>pgd_north:
            pgd_north=max(abs(n[0].data))
        if max(abs(e[0].data))>pgd_east:
            pgd_east=max(abs(e[0].data))
        if max(abs(z[0].data))>pgd_up:
            pgd_up=max(abs(z[0].data))
        
        if k==0:
            N=n.copy()
            E=e.copy()
            Z=z.copy()
        else:
            N+=n
            E+=e
            Z+=z
        
        #Get ray arrivals
        deg=locations2degrees(hypo[1],hypo[0],lon[k],lon[k])
    
    #Reshape
    RGB=zeros((len(sta),N[0].stats.npts,3))
    for k in range(len(sta)):
        
        #normalize and scale
        #n=abs(N[k].data)/pgd_north
        #e=abs(E[k].data)/pgd_east
        #z=abs(Z[k].data)/pgd_up
        
        n=log10(abs(N[k].data))-log10(lower_bound)
        e=log10(abs(E[k].data))-log10(lower_bound)
        z=log10(abs(Z[k].data))-log10(lower_bound)
        n=n/n.max()
        e=e/e.max()
        z=z/z.max()
        

        
        RGB[k,:,0]=abs(n-1)
        RGB[k,:,1]=abs(e-1)
        RGB[k,:,2]=abs(z-1)
        
        #RGB[:,k,0]=(n/2.)+0.5
        #RGB[:,k,1]=(e/2.)+0.5
        #RGB[:,k,2]=(z/2.)+0.5
    
        
 
    plt.figure()
    plt.imshow(RGB,aspect='auto',origin='lower')
    plt.xlim(tlims)
    plt.xlabel('Time (s)')
    plt.ylabel('Station No.')
    plt.title(rupture+', Mw = '+str(round(Mw,2)),fontsize=14)
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
    plt.scatter(sta[iup,1],sta[iup,2],c='r',s=sta[iup,5]*verti,lw=0.1)
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
        print('Reading data from rupture '+ruptures[krupt])
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
    
    
    
def plot_3stations(home,project_name,run_name,run_number,GF_list,stations=['chzz','p154','tfno']):
    '''
    Plot 3 sample waveforms
    '''
    
    from obspy.geodetics.base import gps2dist_azimuth
    from numpy import genfromtxt,array,ones,zeros,where
    from matplotlib import pyplot as plt
    from string import replace
    from obspy import read
    from obspy.taup import TauPyModel
    from obspy.geodetics import locations2degrees
    from matplotlib.ticker import MultipleLocator
    
    xmajorLocator = MultipleLocator(50)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(200)
    yminorLocator = MultipleLocator(50)
    
    # Read summary file
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    event_log=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
    i1=where(sta==stations[0])[0][0]
    i2=where(sta==stations[1])[0][0]
    i3=where(sta==stations[2])[0][0]
    i=[i1,i2,i3]
    sta=sta[i]
    lonlat=lonlat[i]
    
    #Load velocity model for ray tracing
    velmod = TauPyModel(model="/Users/dmelgar/FakeQuakes/Cascadia/structure/cascadia")
    
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
    
    #Get p and s phases
    ptime=1e6*ones(len(sta))
    stime=1e6*ones(len(sta))
    for k in range(len(sta)):        
        #Get ray arrivals
        deg=locations2degrees(hypo[1],hypo[0],lonlat[k,1],lonlat[k,0])
        # Ray trace
        arrivals = velmod.get_travel_times(source_depth_in_km=hypo[2],distance_in_degree=deg,phase_list=['P','Pn','S','Sn','p','s'])
        #Determine P and S arrivals
        for kphase in range(len(arrivals)):
            if 'P' == arrivals[kphase].name or 'p' == arrivals[kphase].name or 'Pn' == arrivals[kphase].name:
                if arrivals[kphase].time<ptime[k]:
                    ptime[k]=arrivals[kphase].time
            if 'S' == arrivals[kphase].name or 's' == arrivals[kphase].name or 'Sn' == arrivals[kphase].name:
                if arrivals[kphase].time<stime[k]:
                    stime[k]=arrivals[kphase].time
                    
    plt.figure(figsize=[16,4])
    
    plt.subplot(131)
    eoffset=3 ; noffset=0 ; zoffset=-1.5
    n=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[0]+'.LYN.sac')
    e=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[0]+'.LYE.sac')
    z=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[0]+'.LYZ.sac')
    plt.plot(e[0].times(),e[0].data+eoffset,'k',lw=1,label='East')
    plt.plot(n[0].times(),n[0].data+noffset,'#FF6347',lw=1,label='North')
    plt.plot(z[0].times(),z[0].data+zoffset,'#4169E1',lw=1,label='Vertical')
    plt.xlabel('Seconds after OT')
    plt.ylabel('Displacement (m)')
    plt.legend(frameon=False)
    plt.xlim([0,200])
    plt.title('Station '+sta[0].upper())
    xmajorLocator = MultipleLocator(50)
    xminorLocator = MultipleLocator(5)
    ymajorLocator = MultipleLocator(1)
    yminorLocator = MultipleLocator(0.2)
    ax=plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)   
    yl=ax.get_ylim()
    plt.plot([ptime[0],ptime[0]],yl,'--',lw=3,c='#A0A0A0')
    plt.plot([stime[0],stime[0]],yl,'--',lw=3,c='#A0A0A0')
    
    plt.subplot(132)
    eoffset=1 ; noffset=0 ; zoffset=-0.5
    n=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[1]+'.LYN.sac')
    e=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[1]+'.LYE.sac')
    z=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[1]+'.LYZ.sac')
    plt.plot(e[0].times(),e[0].data+eoffset,'k',lw=1,label='East')
    plt.plot(n[0].times(),n[0].data+noffset,'#FF6347',lw=1,label='North')
    plt.plot(z[0].times(),z[0].data+zoffset,'#4169E1',lw=1,label='Vertical')
    plt.xlabel('Seconds after OT')
    plt.xlim([0,300])
    plt.title('Station '+sta[1].upper())
    xmajorLocator = MultipleLocator(50)
    xminorLocator = MultipleLocator(5)
    ymajorLocator = MultipleLocator(0.5)
    yminorLocator = MultipleLocator(0.1)
    ax=plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    yl=ax.get_ylim()
    plt.plot([ptime[1],ptime[1]],yl,'--',lw=3,c='#A0A0A0')
    plt.plot([stime[1],stime[1]],yl,'--',lw=3,c='#A0A0A0')  

    plt.subplot(133)
    eoffset=0.2 ; noffset=0 ; zoffset=-0.1
    n=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[2]+'.LYN.sac')
    e=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[2]+'.LYE.sac')
    z=read(home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/'+sta[2]+'.LYZ.sac')
    plt.plot(e[0].times(),e[0].data+eoffset,'k',lw=1,label='East')
    plt.plot(n[0].times(),n[0].data+noffset,'#FF6347',lw=1,label='North')
    plt.plot(z[0].times(),z[0].data+zoffset,'#4169E1',lw=1,label='Vertical')
    plt.xlabel('Seconds after OT')
    plt.xlim([0,500])
    plt.title('Station '+sta[2].upper())
    xmajorLocator = MultipleLocator(100)
    xminorLocator = MultipleLocator(10)
    ymajorLocator = MultipleLocator(0.1)
    yminorLocator = MultipleLocator(0.02)
    ax=plt.gca()
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)  
    yl=ax.get_ylim()
    plt.plot([ptime[2],ptime[2]],yl,'--',lw=3,c='#A0A0A0')
    plt.plot([stime[2],stime[2]],yl,'--',lw=3,c='#A0A0A0')
    
    plt.subplots_adjust(left=0.05,right=0.98,top=0.9,bottom=0.2,wspace=0.17)
    
    
    
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
        print(paths[krupt])
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
   
            
                              
            
def source_time_function(rupt,epicenter,dt=0.001,t_total=500,stf_type='dreger',plot=True):
    '''
    Plot source time function of complete rupture
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,unique,log10,where,floor,argmin,r_,zeros
    from mudpy.forward import build_source_time_function
    
    f=genfromtxt(rupt)
    num=f[:,0]
    #Get slips
    ss=f[:,8]
    ds=f[:,9]
    #Add it up
    slip=(ss**2+ds**2)**0.5
    #Now parse for multiple rupture speeds
    unum=unique(num)
    nfault=len(unum)
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Get rupture times for subfault windows
    trup=f[:,12]
    #Loop over subfaults
    faults_added=0
    for kfault in range(nfault):
#        if kfault%50==0:
#            print('Working on subfault %d of %d' % (kfault,nfault))
        if rise_time[kfault]>0:
            faults_added+=1
            #get stf
            t,Mdot=build_source_time_function(rise_time[kfault],dt,t_total,stf_type=stf_type)
            #What is the moment at this subfault?
            moment=slip[kfault]*mu[kfault]*area[kfault]
            #Scale stf by this moment
            scale_factor=moment/dt
            Mdot=Mdot*scale_factor
            #Shift according to rupture onset time
            tdiff=t-trup[kfault]
            ishift=argmin(abs(tdiff))
            if ishift!=0:
                Mdot=r_[zeros(ishift),Mdot[:-ishift]]
            if faults_added==1:
                Mrate=Mdot.copy()
            else:
                Mrate+=Mdot

    #Get power
    exp=floor(log10(Mrate.max()))
    M1=Mrate/(10**exp)
    if plot==True:
        plt.figure()
        plt.fill(t,Mrate,'b',alpha=0.5)
        plt.plot(t,Mrate,color='k',lw=0.1)
        plt.grid()
        plt.xlabel('Time(s)')
        #plt.ylabel('Moment Rate ('+r'$\times 10^{'+str(int(exp))+r'}$Nm/s)')
        plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.7, wspace=0, hspace=0)
        plt.show()
    return t,Mrate  
    
         
def slip_rate_and_spectrum(rupt,epicenter,dt=0.005,t_total=50,ref1offset=-10,ref2offset=-10):
    '''
    Plot average slip rate function and spectrum
    '''
    import matplotlib.pyplot as plt
    from numpy import genfromtxt,unique,log10,where,floor,argmin,r_,zeros,logspace
    from mudpy.forward import build_source_time_function
    import nitime.algorithms as tsa
    from scipy.signal import periodogram
    from matplotlib.ticker import MultipleLocator
    
    f=genfromtxt(rupt)
    num=f[:,0]
    #Get slips
    ss=f[:,8]
    ds=f[:,9]
    #Add it up
    slip=(ss**2+ds**2)**0.5
    #Now parse for multiple rupture speeds
    unum=unique(num)
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    #Moment
    moment=mu*area*slip
    #Get rise times
    rise_time=f[0:len(unum),7]
    i=where(rise_time>0)[0]
    mean_rise_time=rise_time[i].mean()
    mean_moment=moment[i].mean()
    
    #Get stfs
    t,Mdot1=build_source_time_function(mean_rise_time,dt,t_total,stf_type='triangle')
    t,Mdot2=build_source_time_function(mean_rise_time,dt,t_total,stf_type='dreger',zeta=0.2)
    t,Mdot3=build_source_time_function(mean_rise_time,dt,t_total,stf_type='dreger',zeta=0.1)
    Mdot1=r_[zeros(2000),Mdot1[0:-2000]]
    Mdot2=r_[zeros(2000),Mdot2[0:-2000]]
    
    Mdot1=Mdot1/dt
    Mdot2=Mdot2/dt
    Mdot3=Mdot3/dt
    
    
    #Get psd
    #f, p_triangle, nu = tsa.multi_taper_psd(Mdot1,Fs=1./dt,adaptive=False,jackknife=False,low_bias=True)
    #f, p_dreger, nu = tsa.multi_taper_psd(Mdot2,Fs=1./dt,adaptive=False,jackknife=False,low_bias=True)
    
    fp,pp_triangle=periodogram(Mdot1, fs=1./dt, nfft=None, return_onesided=False,window='hamming',scaling='spectrum')
    fp,pp_dreger=periodogram(Mdot3, fs=1./dt, nfft=None, return_onesided=False,window=None,scaling='spectrum')
    
    
    plt.figure(figsize=(9,4))
    
    plt.subplot(121)
    plt.plot(t,Mdot1,c='#DC143C',lw=2,label='Triangle STF')
    plt.plot(t,Mdot2,c='#606060',lw=2,label='Dreger STF')
    plt.xlim([0,50])
    plt.ylim([0,0.18])
    plt.xlabel('Seconds')
    plt.ylabel('Moment rate')
    ax=plt.gca()
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    plt.legend(frameon=False,ncol=2,labelspacing=0.1,bbox_to_anchor=(1.8, 1.18))
    ax.set_xlim([0,50])
    xmajorLocator = MultipleLocator(10)
    xminorLocator = MultipleLocator(2)
    ymajorLocator = MultipleLocator(0.05)
    yminorLocator = MultipleLocator(0.01)
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)






    plt.subplot(122)
    plt.loglog(fp,pp_triangle**0.5,c='#DC143C')
    plt.loglog(fp,pp_dreger**0.5,c='#606060',lw=2)  
    plt.ylim([5e-8,1e-1])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Fourier amplitude')
    ax=plt.gca()
    
    
    fref=logspace(0,1)
    ref1=10**(-log10(fref)+ref1offset)
    ref2=10**(-2*log10(fref)+ref2offset)
    plt.loglog(fref,ref1,'--',c='#1E90FF',lw=3)
    plt.loglog(fref,ref2,'--',c='#1E90FF',lw=3)
    plt.annotate(r'$f^{-1}$', xy=(3.5, 0.0006))
    plt.annotate(r'$f^{-2}$', xy=(3.5, 0.000006))
    
    
    plt.xlim([0.02,10])
    ax=plt.gca()
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.tick_params(which='major',length=7,width=1)
    ax.tick_params(which='minor',length=4,width=1)
    
    
    plt.subplots_adjust(bottom=0.16,top=0.9,left=0.1,right=0.9,wspace=0.05)
    
    plt.show()


def fence_slip(home,project_name,run_name,run_number,maxslip=None,UTM_zone='10S',elev=20,azimuth=None,fudge=10,fake_hypo=[0.1,0.1],
        borderwidth=0.5,figsize=(21,7),xtick=10,ytick=10,ztick=5,inverse_model=False,hypocenter=None):
    '''
    Make fence diagram of rupture file
    '''
    
    from numpy import genfromtxt,array,zeros,where
    from matplotlib import pyplot as plt
    from matplotlib import colors
    from pyproj import Proj
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from string import replace
    from matplotlib import ticker
    from matplotlib.ticker import MultipleLocator
    from mudpy import gmttools

    #Get rupture data
    if inverse_model==False:
        fault=genfromtxt(home+project_name+'/output/ruptures/%s.%s.rupt' % (run_name,run_number))
        #Parse log file for hypocenter
        log_file=home+project_name+'/output/ruptures/%s.%s.log' % (run_name,run_number)
        f=open(log_file,'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Hypocenter (lon,lat,z[km])' in line:                
                s=replace(line.split(':')[-1],'(','')
                s=replace(s,')','')
                hypocenter=array(s.split(',')).astype('float')
                loop_go=False  
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])   
        f.close() 
    else:
        fault_name=home+project_name+'/output/inverse_models/models/%s.%s.inv' % (run_name,run_number)
        gmttools.make_total_model(fault_name,thresh=0)
        fault=genfromtxt(home+project_name+'/output/inverse_models/models/%s.%s.inv.total' % (run_name,run_number))
        #Parse log file for hypocenter
        log_file=home+project_name+'/output/inverse_models/models/%s.%s.log' % (run_name,run_number)
        f=open(log_file,'r')
        loop_go=True
        while loop_go:
            line=f.readline()  
            if 'Mw' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])   
                loop_go=False
        f.close() 
    


    
    #get subfault corners
    corners=get_subfault_corners(fault)
    
    #Convert ot UTM (in km)
    corners=corners2utm(corners,UTM_zone=UTM_zone)
    
    #Get UTM coords of hypocenter
    P=Proj(proj='utm',ellps='WGS84',zone=UTM_zone)
    hypocenter[0],hypocenter[1]=P(hypocenter[0],hypocenter[1])
    hypocenter[0]/=1000
    hypocenter[1]/=1000
    
    #Make hypocenter the origin
    corners[:,0]-=hypocenter[0]
    corners[:,3]-=hypocenter[0]
    corners[:,6]-=hypocenter[0]
    corners[:,9]-=hypocenter[0]
    corners[:,1]-=hypocenter[1]
    corners[:,4]-=hypocenter[1]
    corners[:,7]-=hypocenter[1]
    corners[:,10]-=hypocenter[1]
    
    #Get mean strike for inital viewing angle
    strike=fault[:,4].mean()
    
    #Normalized slip
    slip=(fault[:,8]**2+fault[:,9]**2)**0.5
    
    #Saturate to maxslip
    if maxslip!=None:
        imax=where(slip>maxslip)[0]
        slip[imax]=maxslip
    #normalize
    norm_slip=slip/slip.max()
    
    #Get colormaps
    pqlx = colors.LinearSegmentedColormap('pqlx',pqlx_dict,256)
    whitejet =colors.LinearSegmentedColormap('whitejet',whitejet_dict,256)
    
    #Azimuth viewing angle
    if azimuth==None:
        azimuth=strike+90
    
    #Plot init, axes positions etc
    fig=plt.figure(figsize=figsize)
    
    ax1 = fig.add_subplot(311, projection='3d')
    ax1.set_xlim([corners[:,0].min()-fudge,corners[:,0].max()+fudge])
    ax1.set_ylim([corners[:,1].min()-fudge,corners[:,1].max()+fudge])
    ax1.set_zlim([corners[:,2].min()-fudge/4,corners[:,2].max()+fudge/4])
    #Fenagle the axis ticks
    xmajorLocator = MultipleLocator(xtick)
    ymajorLocator = MultipleLocator(ytick)
    zmajorLocator = MultipleLocator(ztick)
    ax1.xaxis.set_major_locator(xmajorLocator)
    ax1.yaxis.set_major_locator(ymajorLocator)
    ax1.zaxis.set_major_locator(zmajorLocator)
    ax1.invert_zaxis()
    ax1.view_init(elev=elev, azim=azimuth)

    #Make one patch per subfault
    for ksub in range(len(corners)):
        vertices=[[tuple(corners[ksub,0:3]),tuple(corners[ksub,3:6]),tuple(corners[ksub,6:9]),tuple(corners[ksub,9:12])]]
        subfault=Poly3DCollection(vertices, linewidths=borderwidth)
        #subfault.set_color(pqlx(norm_slip[ksub]))
        subfault.set_color(pqlx(norm_slip[ksub]))
        subfault.set_linewidth(borderwidth)
        subfault.set_edgecolor('#505050')
        ax1.add_collection3d(subfault)
    
    #Hypocenter
    ax1.scatter(fake_hypo[0],fake_hypo[1],hypocenter[2],s=220,marker='*',c='#7CFC00')
    
    #Dummy mapable for colorbar
    s=plt.scatter(zeros(len(fault)),zeros(len(fault)),c=slip,cmap=pqlx,s=0.00001,lw=0)
    
    #Mke colorbar
    cb=plt.colorbar(s,shrink=0.9,pad=-0.07)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    cb.set_label('Slip (m)')
    
    #Labels n' stuff
    ax1.set_xlabel('\n\nEast (km)')
    ax1.set_ylabel('\n\nNorth (km)')
    ax1.set_zlabel('Depth (km)',rotation=90)
    #plt.title(home+project_name+'/output/ruptures/%s.%s.rupt' % (run_name,run_number))
    plt.title(run_name+' '+run_number+' Mw '+str(Mw))

    




    # Ok now let's do the rise times
    
    rise=fault[:,7]
    norm_rise=rise/rise.max()
    
    #Get colormaps
    risemap = colors.LinearSegmentedColormap('risetimes',rise_times_dict,256)
    
    #Azimuth viewing angle
    if azimuth==None:
        azimuth=strike+90
    
    #Plot init, axes positions etc
    ax2 = fig.add_subplot(312, projection='3d')
    ax2.set_xlim([corners[:,0].min()-fudge,corners[:,0].max()+fudge])
    ax2.set_ylim([corners[:,1].min()-fudge,corners[:,1].max()+fudge])
    ax2.set_zlim([corners[:,2].min()-fudge/4,corners[:,2].max()+fudge/4])
    #Fenagle the axis ticks
    ax2.xaxis.set_major_locator(xmajorLocator)
    ax2.yaxis.set_major_locator(ymajorLocator)
    ax2.zaxis.set_major_locator(zmajorLocator)
    ax2.invert_zaxis()
    ax2.view_init(elev=elev, azim=azimuth)

    #Make one patch per subfault
    for ksub in range(len(corners)):
        vertices=[[tuple(corners[ksub,0:3]),tuple(corners[ksub,3:6]),tuple(corners[ksub,6:9]),tuple(corners[ksub,9:12])]]
        subfault=Poly3DCollection(vertices, linewidths=borderwidth)
        subfault.set_color(risemap(norm_slip[ksub]))
        subfault.set_linewidth(borderwidth)
        subfault.set_edgecolor('#505050')
        ax2.add_collection3d(subfault)

    #Hypocenter
    ax2.scatter(fake_hypo[0],fake_hypo[1],hypocenter[2],s=220,marker='*',c='#7CFC00')    
            
    #Dummy mapable for colorbar
    s=plt.scatter(zeros(len(fault)),zeros(len(fault)),c=rise,cmap=risemap,s=0.00001,lw=0)
    
    #Mke colorbar
    cb=plt.colorbar(s,shrink=0.9,pad=-0.07)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    cb.set_label('Rise time (s)')
    
    #Labels n' stuff
    ax2.set_xlabel('\n\nEast (km)')
    ax2.set_ylabel('\n\nNorth (km)')
    ax2.set_zlabel('Depth (km)',rotation=90)  




    # Finally the onset times
    
    onset=fault[:,12]
    norm_onset=onset/onset.max()
    
    #Get colormaps
    onsetmap = colors.LinearSegmentedColormap('whitejet',whitejet_dict,256)
    
    #Azimuth viewing angle
    if azimuth==None:
        azimuth=strike+90
    
    #Plot init, axes positions etc
    ax3 = fig.add_subplot(313, projection='3d')
    ax3.set_xlim([corners[:,0].min()-fudge,corners[:,0].max()+fudge])
    ax3.set_ylim([corners[:,1].min()-fudge,corners[:,1].max()+fudge])
    ax3.set_zlim([corners[:,2].min()-fudge/4,corners[:,2].max()+fudge/4])
    #Fenagle the axis ticks
    ax3.xaxis.set_major_locator(xmajorLocator)
    ax3.yaxis.set_major_locator(ymajorLocator)
    ax3.zaxis.set_major_locator(zmajorLocator)
    ax3.invert_zaxis()
    ax3.view_init(elev=elev, azim=azimuth)

    #Make one patch per subfault
    for ksub in range(len(corners)):
        vertices=[[tuple(corners[ksub,0:3]),tuple(corners[ksub,3:6]),tuple(corners[ksub,6:9]),tuple(corners[ksub,9:12])]]
        subfault=Poly3DCollection(vertices, linewidths=borderwidth)
        subfault.set_color(onsetmap(norm_onset[ksub]))
        #subfault.set_color(plt.cm.jet(norm_onset[ksub]))
        subfault.set_linewidth(borderwidth)
        subfault.set_edgecolor('#505050')
        ax3.add_collection3d(subfault)

    #Hypocenter
    ax3.scatter(fake_hypo[0],fake_hypo[1],hypocenter[2],s=220,marker='*',c='#7CFC00')  
         
    #Dummy mapable for colorbar
    s=plt.scatter(zeros(len(fault)),zeros(len(fault)),c=onset,cmap=onsetmap,s=0.00001,lw=0)
    #s=plt.scatter(zeros(len(fault)),zeros(len(fault)),c=onset,cmap=plt.cm.jet,s=0.00001,lw=0)
    
    #Mke colorbar
    cb=plt.colorbar(s,shrink=0.9,pad=-0.07)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cb.locator=tick_locator
    cb.update_ticks()
    cb.set_label('Onset time (s)')
    
    #Labels n' stuff
    ax3.set_xlabel('\n\nEast (km)')
    ax3.set_ylabel('\n\nNorth (km)')
    ax3.set_zlabel('Depth (km)',rotation=90)          
      
    
                  
    #Fix the spacing in between subplots
    plt.subplots_adjust(left=0.01,right=0.99,bottom=0.03,top=0.98,hspace=0.05)               
                                                         
    plt.show()

    
    

    

def get_subfault_corners(fault):
    '''
    Calculate the coordiantes of the 4 corners making up every subfault
    '''
    
    from pyproj import Geod
    from numpy import cos,deg2rad,sin,zeros
    
    #Projection object
    g=Geod(ellps='WGS84')
    
    #Output variable
    corners=zeros((len(fault),12))
    
    for kfault in range(len(fault)):
        lon=fault[kfault,1]
        lat=fault[kfault,2]
        z=fault[kfault,3]
        W=fault[kfault,11]
        L=fault[kfault,10]
        strike=fault[kfault,4]
        dip=fault[kfault,5]
        
        #Get top center and bottom center
        delta_h=W*cos(deg2rad(dip))/2
        lon_top_center,lat_top_center,foo=g.fwd(lon,lat,strike-90,delta_h)
        lon_bot_center,lat_bot_center,foo=g.fwd(lon,lat,strike+90,delta_h)
        z_top=z-(W*sin(deg2rad(dip))/2)/1000
        z_bot=z+(W*sin(deg2rad(dip))/2)/1000
        
        #Get coordinates of corners
        
        #Corner1
        lon_corner,lat_corner,foo=g.fwd(lon_top_center,lat_top_center,strike,L/2)
        corners[kfault,0:3]=lon_corner,lat_corner,z_top
        
        #Corner2
        lon_corner,lat_corner,foo=g.fwd(lon_top_center,lat_top_center,strike+180,L/2)
        corners[kfault,3:6]=lon_corner,lat_corner,z_top
        
        #Corner3
        lon_corner,lat_corner,foo=g.fwd(lon_bot_center,lat_bot_center,strike+180,L/2)
        corners[kfault,6:9]=lon_corner,lat_corner,z_bot
        
        #Corner4
        lon_corner,lat_corner,foo=g.fwd(lon_bot_center,lat_bot_center,strike,L/2)
        corners[kfault,9:12]=lon_corner,lat_corner,z_bot
        
    return corners
    
def corners2utm(corners,UTM_zone='10S'):
    '''
    Convert coordiantes of the corners to UTM
    '''
    
    from pyproj import Proj

    P=Proj(ellps='WGS84',proj='utm',zone=UTM_zone)
    corners_out=corners.copy()
    
    corners_out[:,0],corners_out[:,1]=P(corners[:,0],corners[:,1])
    corners_out[:,3],corners_out[:,4]=P(corners[:,3],corners[:,4])
    corners_out[:,6],corners_out[:,7]=P(corners[:,6],corners[:,7])
    corners_out[:,9],corners_out[:,10]=P(corners[:,9],corners[:,10])
    
    corners_out[:,0]/=1000
    corners_out[:,1]/=1000
    corners_out[:,3]/=1000
    corners_out[:,4]/=1000
    corners_out[:,6]/=1000
    corners_out[:,7]/=1000
    corners_out[:,9]/=1000
    corners_out[:,10]/=1000
    
    return corners_out
    
    
def plot_dtopo(dtopo_file,s=5):
    
    from matplotlib import pyplot as plt
    from numpy import where,genfromtxt
    
    d=genfromtxt(dtopo_file)
    i=where(d[:,0]==d[:,0].max())[0]
    d=d[i,:]
    
    max_vert=max(abs(d[:,3]))
    
    plt.figure()
    plt.scatter(360+d[:,1],d[:,2],c=d[:,3],vmin=-max_vert,vmax=max_vert,cmap=plt.cm.seismic)
    plt.colorbar(label='Vertical deformation (m)')
    plt.title(dtopo_file)
    plt.show()
    
    
def plot_hypocenter_locations(home,project_name,run_name):
    
    from matplotlib import pyplot as plt
    from glob import glob
    from numpy import array
    
    logs=glob(home+project_name+'/output/ruptures/'+run_name+'*.log')
    hypo_lon=[]
    hypo_lat=[]
    hypo_z=[]
    centroid_lon=[]
    centroid_lat=[]
    centroid_z=[]
    L=[]
    for kfault in range(len(logs)):
        
        #Get info about fault
        f=open(logs[kfault],'r')
        loop_go=True
        while loop_go:
            line=f.readline()
            if 'Maximum length Lmax:' in line:
                l=float(line.split(':')[-1].split()[0])
                L.append(l)
            if 'Hypocenter (lon,lat,z[km])' in line:                
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
                hypo=array(s.split(',')).astype('float')
                
                hypo_lon.append(hypo[0])
                hypo_lat.append(hypo[1])
                hypo_z.append(hypo[2])
                
                
            if 'Centroid (lon,lat,z[km])' in line:                
                s=line.split(':')[-1]
                s=s.replace('(','')
                s=s.replace(')','')
                centroid=array(s.split(',')).astype('float')
                
                centroid_lon.append(centroid[0])
                centroid_lat.append(centroid[1])
                centroid_z.append(centroid[2])
     
                loop_go=False  
                
#    plt.figure()
#    plt.scatter(hypo_lon,hypo_lat,c=hypo_z,cmap='jet')
#    plt.colorbar(label='Hypocenter Depth (km)')
#    plt.axis('equal')
#    
#    plt.figure()
#    plt.scatter(centroid_lon,centroid_lat,c=hypo_z,cmap='jet')
#    plt.colorbar(label='Centroid Depth (km)')
#    plt.axis('equal')
    
    plt.figure()
    plt.scatter(hypo_lon,hypo_lat,facecolor='k')
    plt.scatter(centroid_lon,centroid_lat,facecolor='r')
    plt.legend(['Hypocenters','Centroids'])
    for k in range(len(hypo_lon)):
        plt.plot([hypo_lon[k],centroid_lon[k]],[hypo_lat[k],centroid_lat[k]],'b')
    plt.axis('equal')
        
    plt.figure()
    for k in range(len(hypo_lon)):
        plt.scatter(hypo_lat[k],hypo_lat[k]-centroid_lat[k])
    plt.ylabel('Hypo lat - centroid lat')
    
    plt.figure()
    for k in range(len(hypo_lon)):
        plt.scatter(L[k],abs(hypo_lat[k]-centroid_lat[k])*110)
    plt.ylabel('Hypo lat - centroid lat')
    plt.xlabel('Length (km)')
    
    plt.show()
