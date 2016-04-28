'''
Module for plotting FakeQuakes related studd
'''
from matplotlib import rc
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

rc('font', **font)


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
    #summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
    summary_file=summary_file=home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'

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
    ax.scatter(d,pgd,s=60, facecolors='none', edgecolors='m',marker='d',label='cascadia.mod')
    ax.scatter(d,pgd2,s=60, facecolors='none', edgecolors='k',marker='d',label='gil7.mod')
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
        summary_file=summary_file=home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
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
        summary_file='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
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
        for k in range(len(bin_edges)-1):
            i=where((d_all>=bin_edges[k]) & (d_all<bin_edges[k+1]))[0]
            misfit_bin[k]=misfit[i].mean()
            misfit_bin2[k]=misfit2[i].mean()
            bin_centers[k]=10**(log10(bin_edges[k+1])-(log10(bin_edges[k+1])-log10(bin_edges[k]))/2)
        
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(bin_centers,misfit_bin,c='k',lw=2,label='cascadia.mod')
        ax.plot(bin_centers,misfit_bin2,c='r',lw=2,label='gil7.mod')
        plt.legend()
        ax.scatter(bin_centers,misfit_bin,c='k',lw=0.5,s=70)
        ax.scatter(bin_centers,misfit_bin2,c='r',lw=0.5,s=70)
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
                                                        
 
 
def pgd_magnitude_misfit(home,project_name,hist='True',rupture_list='ruptures.list',Mw_lims=[7.8,9.5],A=-4.434,B=1.047,C=-0.138,misfit_lims=[-3,3],nbins=10):
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
        summary_file=summary_file=home+project_name+'/output/cosine/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
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
        summary_file='/Users/dmelgar/FakeQuakes/Cascadia_gil7/output/cascadia.mod_scenarios/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
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
        plt.xlim(Mw_lims)
    
        plt.ylim(misfit_lims)
        plt.subplots_adjust(bottom=0.15)
    else:
        bin_edges=linspace(Mw_all.min(),Mw_all.max(),nbins)
        bin_centers=zeros(len(bin_edges)-1)
        misfit_bin=zeros(len(bin_edges)-1)
        misfit_bin2=zeros(len(bin_edges)-1)
        for k in range(len(bin_edges)-1):
            i=where((Mw_all>=bin_edges[k]) & (Mw_all<bin_edges[k+1]))[0]
            misfit_bin[k]=misfit[i].mean()
            misfit_bin2[k]=misfit2[i].mean()
            bin_centers[k]=bin_edges[k+1]-((bin_edges[k+1]-bin_edges[k])/2)
        
        fig = plt.figure()
        ax = plt.gca()
        ax.plot(bin_centers,misfit_bin,c='k',lw=2,label='cascadia.mod')
        ax.plot(bin_centers,misfit_bin2,c='r',lw=2,label='gil7.mod')
        plt.legend()
        ax.scatter(bin_centers,misfit_bin,c='k',lw=0.5,s=70)
        ax.scatter(bin_centers,misfit_bin2,c='r',lw=0.5,s=70)
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
    
