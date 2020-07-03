# -*- coding: utf-8 -*-
'''
C.Ruhl 08/2016

Functions to support Fakequakes by building fault models, generating response spectra for validation, etc.
'''
import numpy as np
from matplotlib import cm

def bulk_2Dmisfit(home,project_names,rupture_list='ruptures.list',val='SA',Mw_lims=[5.8,9.5],dist_lims=[10,1000],cmapwf=cm.magma_r,misfit_lims=[-3,3],numwf_lims=[0,1000],GOF_lims=[0,2],n_mag_bins=10,n_dist_bins=10,A=-4.434,B=1.047,C=-0.138):
    '''
    Plot spectral acceleration misfit as a function of both distance and magnitude
    val = 'SA' --> Spectral Acceleration (10s period)
    val = 'PGD' --> Peak Ground Displacement 
    '''
    
    import numpy as np
    from numpy import genfromtxt,array,zeros,logspace,linspace,r_,log,genfromtxt,log10,ones,where,arange,mean
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import MaxNLocator
    from string import replace
    #import pylab
    from obspy.geodetics.base import gps2dist_azimuth
    
    if val=='SA':
        SAobs_all=array([])
        rjb_all=array([])
        Mw_all=array([])
        SAcalc_all=array([])
    elif val=='PGD':
        pgd_all=array([])
        d_all=array([])
        Mw_all=array([])
        pgd_predicted_all=array([])
    
    for i in range(len(project_names)):
        project_name=project_names[i]
        
        if project_name=='Cascadia':
            ruptures=genfromtxt(home+project_name+'/'+rupture_list,dtype='S')
        else:
            ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
        for k in range(len(ruptures)):
            run_name=ruptures[k].split('.')[0]
            run_number=ruptures[k].split('.')[1]
            
            if val=='SA':    
                # Read analysis file
                if project_name=='Cascadia':
                    analysis_file=home+project_name+'/'+run_name+'.'+run_number+'/_analysis.'+run_name+'.'+run_number+'.txt'
                else:
                    analysis_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_analysis.'+run_name+'.'+run_number+'.txt'
                SAobs=genfromtxt(analysis_file,usecols=[5])      
                SAcalc=genfromtxt(analysis_file,usecols=[6])
                rjb=genfromtxt(analysis_file,usecols=[8])            
                # get the magnitude of each rupture from the log file
                if project_name=='Cascadia':
                    logfile=home+project_name+'/'+run_name+'.'+run_number+'/_'+run_name+'.'+run_number+'.log'
                else:
                    logfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
                f=open(logfile,'r')
                loop_go=True
                while loop_go:
                    line=f.readline()  
                    if 'Actual magnitude' in line:
                        Mw=float(line.split(':')[-1].split(' ')[-1])
                        break
                
                #Concatente to output variables
                SAobs_all=r_[SAobs_all,SAobs]
                rjb_all=r_[rjb_all,rjb]
                Mw_all=r_[Mw_all,Mw*ones(len(rjb))]                   
                SAcalc_all=r_[SAcalc_all,SAcalc]
                
            elif val=='PGD':
                # Read summary file
                if project_name=='Cascadia':
                    summary_file=summary_file=home+project_name+'/'+run_name+'.'+run_number+'/_'+run_name+'.'+run_number+'.offsets'
                else:
                    summary_file=summary_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
                    
                # summary_file=summary_file=home+project_name+'/output/dreger_4tau/'+run_name+'.'+run_number+'/_summary.'+run_name+'.'+run_number+'.txt'
                lonlat=genfromtxt(summary_file,usecols=[1,2])
                pgd=genfromtxt(summary_file,usecols=[6])*100
                # Get hypocenter or centroid
                if project_name=='Cascadia':
                    logfile=home+project_name+'/'+run_name+'.'+run_number+'/_'+run_name+'.'+run_number+'.log'
                else:
                    logfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
                f=open(logfile,'r')
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
    if val=='SA':            
        misfit=log(SAobs_all/SAcalc_all)
    else:  
        misfit=log(pgd_all/pgd_predicted_all)              

    #remove extrneous magnitudes
    i=where((Mw_all>=Mw_lims[0]) & (Mw_all<=Mw_lims[1]))[0]
    Mw_all=Mw_all[i]
    misfit=misfit[i]
    if val=='SA':
        rjb_all=rjb_all[i]
    elif val=='PGD':
        rjb_all=d_all[i]

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
            i=where((Mw_all>=bin_edges_x[kx]) & (Mw_all<bin_edges_x[kx+1]) & (rjb_all>=bin_edges_y[ky]) & (rjb_all<bin_edges_y[ky+1]))[0]
            GOF_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit[i]))
            misfit_bin[kx,ky]=misfit[i].mean()
            frequency_bin[kx,ky]=len(misfit[i])
            bin_centers_x[kx]=bin_edges_x[kx+1]-((bin_edges_x[kx+1]-bin_edges_x[kx])/2)        
            bin_centers_y[ky]=bin_edges_y[ky+1]-((bin_edges_y[ky+1]-bin_edges_y[ky])/2) 
            
    plt.figure(figsize=(15,5.5))    
    plt.subplot(131)                           
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,misfit_bin.T,vmin=-1.5,vmax=1.5,cmap=cm.RdBu_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,GOF_bin.T,colors='#505050',levels=[0.7])
    plt.clabel(CS, inline=1, fontsize=14,fmt='%1.1f')
    ax = plt.gca()    
    ax.set_yscale('log')
    plt.xlim([bin_edges_x.min(),bin_edges_x.max()])
    plt.ylim([bin_edges_y.min(),bin_edges_y.max()])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top") 
    plt.xlabel('Magnitude')
    if val=='SA':
        plt.ylabel('Rjb Distance (km)')
    elif val=='PGD':
        plt.ylabel('Distance (km)')
    levs=ax.get_xticks()
    ax.set_xticklabels(levs,rotation=-55)
    #plt.annotate('(a)',xy=(7.9,680),fontsize=16)
    bbox_props = dict(boxstyle="round", fc="w")
    ax.text(7.9, 350, "(a)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    if val=='SA':
        cb.set_label('Ln(10s SA Obs/Calc)')
    elif val=='PGD':
        cb.set_label('Ln(PGD Obs/Calc)')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplot(132)                           
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,GOF_bin.T,vmin=GOF_lims[0],vmax=GOF_lims[1],cmap=cm.afmhot_r)
    CS = plt.contour(bin_centers_x,bin_centers_y,GOF_bin.T,colors='#808080',levels=[0.7])
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
    PM=plt.pcolormesh(bin_edges_x,bin_edges_y,frequency_bin.T,cmap=cmapwf,vmin=numwf_lims[0],vmax=numwf_lims[1]) 
    CS = plt.contour(bin_centers_x,bin_centers_y,GOF_bin.T,colors='#808080',levels=[0.7])
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
    ax.text(7.9, 350, "(c)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('# of waveforms')
    levs=cb.ax.get_xticks()
    cb.ax.set_xticklabels(levs,rotation=-55)
    tick_locator = MaxNLocator(nbins=6)
    cb.locator = tick_locator
    cb.update_ticks()
    
    plt.subplots_adjust(left=0.07,top=0.87,right=0.98,bottom=0.1,wspace=0.06) 
    if val=='SA':
        plt.savefig(home+'bulk_SA_2D_misfit_and_GOF.pdf')
        return SAobs_all,rjb_all,SAcalc_all,Mw_all
    elif val=='PGD':
        plt.savefig(home+'bulk_PGD_2D_misfit_and_GOF.pdf') 
        d_all=rjb_all
        return pgd_all,d_all,pgd_predicted_all,Mw_all      
    
    
    

def plot_resp_spec(home,project_name,run_name,GF_list,fault_name,rupt,frequency_vector,stalist,slip_type='SS',plot_SA=True,Vs30=760):
    '''
    Function to plot 4, 8, 16, or 30 individual response spectral acceleration
    from station list (stalist) with or without theoretical spectral acceleration 
    (plot_SA) over range of frequencies (frequency_vector).
    
    Christine J. Ruhl, September 2016
    '''
    from numpy import genfromtxt
    from obspy import read
    from matplotlib import pyplot as plt
    from matplotlib import mlab
    
    # Read summary file
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    slon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=1,dtype='float')
    slat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=2,dtype='float')
    
    # get station indices
    staind=[]
    for i in range(len(stalist)):
        if stalist[i] in sta:
            ind=mlab.find(sta==stalist[i])
            staind.append(ind[0])
            
    if len(staind) <= 4:
        subR=2
        subC=2
        subcountR=[0,0,1,1]
        subcountC=[0,1,0,1]
    elif len(staind) <=8:
        subR=4
        subC=2
        subcountR=[0,0,1,1,2,2,3,3]
        subcountC=[0,1,0,1,0,1,0,1]
    elif len(staind) <=16:
        subR=4
        subC=4
        subcountR=[0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]
        subcountC=[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
    elif len(staind) <=30:
        subR=6
        subC=5
        subcountR=[0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5]
        subcountC=[0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4]
    else:
        print('ERROR: Station list too long, shorten to less than or equal to 30!')
    
    fig, axarr = plt.subplots(subR,subC,sharex=False,sharey=False)
    plt.suptitle(run_name+' '+rupt)
    
    #for k in range(len(sta)):
    count=-1
    for k in staind:        
        count=count+1
        ax1 = axarr[subcountR[count],subcountC[count]]
        e=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[k]+'.LYE.sac')
        n=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[k]+'.LYN.sac')
        z=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[k]+'.LYZ.sac')
        # diff twice to get acceleration
        e.differentiate()
        e.differentiate()
        n.differentiate()
        n.differentiate()
        z.differentiate()
        z.differentiate()
        # calculate response spectra using acceleration (converted to g)
        respe=responseSpectrum(e[0].stats.sampling_rate,e[0].data/9.81,frequency_vector)
        respn=responseSpectrum(n[0].stats.sampling_rate,n[0].data/9.81,frequency_vector)
        respz=responseSpectrum(z[0].stats.sampling_rate,z[0].data/9.81,frequency_vector)
        ax1.loglog(frequency_vector,respe)
        ax1.loglog(frequency_vector,respn,'--') 
        ax1.loglog(frequency_vector,respz,':')
        # RtD50 response spectrum
        resp=rotatedResponseSpectrum(n[0].stats.sampling_rate, n[0].data/9.81, e[0].data/9.81, frequency_vector)
        ax1.loglog(frequency_vector,resp[0],'-.')
        
        if plot_SA==True:
            #Get info about fault and rupture
            logfile=home+project_name+'/output/ruptures/'+run_name+'.'+rupt+'.log'
            f=open(logfile,'r')
            loop_go=True
            while loop_go:
                line=f.readline()  
                if 'Actual magnitude' in line:
                    Mw=float(line.split(':')[-1].split(' ')[-1])
                    break
            
            # calculate Rjb
            ruptfile=home+project_name+'/output/ruptures/'+run_name+'.'+rupt+'.rupt'
            Rjb=calc_rjb(slon[k],slat[k],ruptfile)   
            ax1.set_title(sta[k]+' Rjb= '+str(round(Rjb))+' km')
            
            # set Vs30
            Vs30=np.array(Vs30) # keep constant?
            
            # calculate SA based on slip_type
            if slip_type=='SS':
                T,sa,std=bssa14(Mw,Rjb,Vs30,U=0,RS=0,NS=0,SS=1,Z1=None,intensity_measure='SA')
            elif slip_type=='NS':
                T,sa,std=bssa14(Mw,Rjb,Vs30,U=0,RS=0,NS=1,SS=0,Z1=None,intensity_measure='SA')
            elif slip_type=='RS':
                T,sa,std=bssa14(Mw,Rjb,Vs30,U=0,RS=1,NS=0,SS=0,Z1=None,intensity_measure='SA')
                       
            #xlim1=ax1.get_xlim()
            #ylim1=ax1.get_ylim()
            ax1.loglog(1/np.array(T),np.array(sa),'k')
            ax1.loglog(1/np.array(T),(np.exp(np.log(sa)-np.log(std))),'k--')
            ax1.loglog(1/np.array(T),(np.exp(np.log(sa)+np.log(std))),'k--')
            
            #ax1.set_xlim(xlim1)
            #ax1.set_ylim(ylim1)
        
        if subcountR[count]==(subR-1):
            ax1.set_xlabel('Frequency (Hz)')        
        if count==0 and plot_SA==True:
            ax1.legend(['e','n','z','rotD50','calc'],loc='lower right')
        elif count==0:
            ax1.legend(['e','n','z','rotD50'],loc='lower right')
        if subcountC[count] == 0:
            ax1.set_ylabel('SA (g)')
               
    plt.show()
    
def calc_rjb(slon,slat,rupt):
    '''
    Short function to calculate nearest epicentral distance to fault trace 
    (with rupture on it) for one station and the rupture file
    
    Gets active subfaults using risetime != 0.
    These are centers of subfaults, so not technically exact Rjb, but should be 
    close enough since subfaults are generally small.    
    
    Christine J. Ruhl, August 2016
    '''
    from numpy import genfromtxt,array
    from pyproj import Geod
    from matplotlib import mlab
    
    # set projection
    g=Geod(ellps='WGS84')    
    # get fault patches that actually have rupture on them using rise times != 0 in column 8 in ruptfile    
    flon=rupt[:,1]
    flat=rupt[:,2]
    frisetime=rupt[:,10]
    
    ind=mlab.find(frisetime)
    flon=flon[ind]
    flat=flat[ind]
    frisetime=frisetime[ind]

    Rjb=[]
    for ksta in range(len(slon)):        
        d=[]
        for ksource in range(len(flon)):
            az,baz,dist=g.inv(flon[ksource],flat[ksource],slon[ksta],slat[ksta])
            d.append(dist/1000)
        Rjb.append(min(d))
    Rjb=array(Rjb)
    return Rjb


import math
from pyproj import Geod
from mudpy import forward

def build_simple_fault(home,faultname,centroid,strike,dip,length,width):
    '''
    Build a simple fault geometry from fault centroid [lon,lat,depth], 
    strike (degrees), dip (degrees), length (km) and width (km).
    Automatically calculates subfault sizes based on full length of fault
    to maintain ~1000 subfaults total.
    
    Christine J. Ruhl, August 2016
    '''
    g=Geod(ellps='WGS84')
    VR=0.8*3600 # 0.8 times Vs ~ 4 km/s
    
    fault=faultname    
    print('working on the '+fault+' fault.')

    midlon=centroid[0]
    midlat=centroid[1]
            
    # build subfault sizes based on total length of fault
    if length<100:
        strikefactor=1
        dipfactor=1
    elif length>=100 and length<200:
        strikefactor=2
        dipfactor=1
    elif length>=200 and length<300:
        strikefactor=2
        dipfactor=2
    elif length>=300 and length<500:
        strikefactor=3
        dipfactor=2
    else:
        strikefactor=4
        dipfactor=4
            
    # build subfaults parameters using strikefactor and dipfactor
    nstrike=int(round((length)/strikefactor))
    if (nstrike % 2 == 0):
        nstrike=nstrike+1
    dx_strike=(length)/nstrike # in km
    num_updip=0
    num_downdip=int(round(width/dipfactor))
    dx_dip=width/num_downdip
    print(dx_dip)
    hypo=np.array([midlon,midlat,dx_dip/2.])            
    risetime=math.sqrt(dx_strike*dx_strike+dx_dip*dx_dip)/VR
        
    outfile=home+'data/model_info/'+fault+'.fault'
    forward.makefault(outfile,strike,dip,nstrike,dx_dip,dx_strike,hypo,num_updip,num_downdip,risetime)

#def build_faults_xml(faults,directories,workingdir,ucerf_xmlfile):
#    '''
#    Build fault(s) from UCERF3 fault_sections.xml file
#    
#    Input:
#    faults is list of fault names direct from UCERF3 (e.g., ['San Jacinto'])
#    directories is list of corresponding directories to put fault info (e.g., ['san_jacinto'])
#    workingdir is full path to where the fault directories will be created
#    ucerf_xmlfile is full path to XML file
#    
#    Notes:
#    This function does not account for overlapping segments.
#    
#    Christine J. Ruhl, August 2016
#    '''
#    import os
#    import math
#    import numpy as np
#    from pyproj import Geod
#    from mudpy import forward
#    import xml.etree.cElementTree as et 
#    from shutil import copyfile as cp
#    
#    g=Geod(ellps='WGS84')
#    VR=0.8*3600 # 0.8 times Vs ~ 4 km/s
#    
#    tree = et.ElementTree(file=ucerf_xmlfile)
#    root=tree.getroot()
#    for i in range(len(faults)):
#        print('working on the '+faults[i]+' fault...')
#        home=workingdir+directories[i]+'/'
#        if not os.path.exists(home):
#            os.mkdir(home)
#        for file in os.listdir(home):
#            if file.endswith(".file"):
#                os.remove(os.path.join(home,file))
#        # reinitialize some stuff        
#        dd=[]
#        dip=[]
#        sd=[] 
#        lon=[]
#        lat=[]
#        newlon=[]
#        newlat=[]
#        tempdist=[]
#        totallength=0
#        secids=[]
#        seccount=0
#        # get list of section IDs first
#        for elem in tree.iterfind('FaultSectionPrefDataList/'):
#            if faults[i] in elem.get('sectionName') and 'instance 0' in elem.get('sectionName'):      
#                if not secids:
#                    secids.append([elem.tag])
#                    dd.append(float(elem.get('dipDirection')))
#                    dip.append(float(elem.get('aveDip')))
#                    sd.append(float(elem.get('aveLowerDepth')))
#                    lat.append([])
#                    lon.append([])
#                    newlon.append([])
#                    newlat.append([])
#                    tempdist.append([])
#                    continue
#                if dd[len(dd)-1]==float(elem.get('dipDirection')) and sd[len(sd)-1]==float(elem.get('aveLowerDepth')) and dip[len(dip)-1]==float(elem.get('aveDip')):
#                    secids[seccount].append(elem.tag)
#                else:
#                    seccount=seccount+1
#                    secids.append([elem.tag])
#                    dd.append(float(elem.get('dipDirection')))
#                    print(elem.get('dipDirection'))
#                    dip.append(float(elem.get('aveDip')))
#                    sd.append(float(elem.get('aveLowerDepth')))
#                    lat.append([])
#                    lon.append([])
#                    newlon.append([])
#                    newlat.append([])
#                    tempdist.append([])
#                                
#        # then go through each section ID to get other information
#        for j in range(len(secids)): # extra loop added to keep track of multiple sections separately     
#            for k in range(len(secids[j])):
#                for elem0 in tree.iterfind('FaultSectionPrefDataList/'+secids[j][k]+'/FaultTrace'):
#               	    for elem1 in elem0:
#           	        # append if it is: 1) First instance, 2) not equal to previous value
#           	        if not lat[j]: #1	    
#                        lat[j].append(float(elem1.get('Latitude')))
#                        lon[j].append(float(elem1.get('Longitude')))
#                        tempdist[j].append(np.nan)
#           	        elif float(elem1.get('Latitude'))!=lat[j][len(lat[j])-1]: #2	    
#                        lat[j].append(float(elem1.get('Latitude')))
#                        lon[j].append(float(elem1.get('Longitude')))                        
#                        az,baz,dist=g.inv(lon[j][len(lon[j])-2],lat[j][len(lat[j])-2],lon[j][len(lon[j])-1],lat[j][len(lat[j])-1],radians=False)
#                        tempdist[j].append(dist/1000)
#                        totallength=totallength+dist/1000
#                            
#        print('length = ',totallength)
#        # calculate Mmax from Blaser et al 2010 relationship
#        # LOG10 (Length) = a + b x Mw
#        a=-2.69
#        b=0.64
#        Mmax=(np.log10(totallength)-a)/b
#        if Mmax>8.5:
#            Mmax=8.5
#        print('Mmax =',Mmax)
#                                            
#        # build subfault sizes based on total length of fault
#        if totallength<100:
#            strikefactor=1
#            dipfactor=1
#        elif totallength>=100 and totallength<150:
#            strikefactor=2
#            dipfactor=1
#        elif totallength>=150 and totallength<250:
#            strikefactor=2
#            dipfactor=2
#        elif totallength>=250 and totallength<400:
#            strikefactor=round(totallength/150)
#            dipfactor=2
#        else:	
#            strikefactor=round(totallength/200)
#            dipfactor=4
#            
#        # loop through distances and get rid of ones less than the size of strikefactor
#        newdist=[]
#        for j in range(len(lat)):
#            for k in range(len(lat[j])):
#                if k==0:
#                    newlat[j].append(lat[j][k])
#                    newlon[j].append(lon[j][k])
#                elif tempdist[j][k]>=(0.75*strikefactor):
#                    newlat[j].append(lat[j][k])
#                    newlon[j].append(lon[j][k])                                        
#                    az,baz,dist=g.inv(lon[j][len(lon[j])-2],lat[j][len(lat[j])-2],lon[j][len(lon[j])-1],lat[j][len(lat[j])-1],radians=False)
#                    newdist.append(dist/1000)
#                    
#        lon=newlon
#        lat=newlat       
#        
#        # check dip direction based on strike 
#        for j in range(len(dd)):                                       
#            az,baz,dist=g.inv(lon[j][len(lon[j])-2],lat[j][len(lat[j])-2],lon[j][len(lon[j])-1],lat[j][len(lat[j])-1],radians=False)
#            if (az<=0):
#                az=360+az
#            assdd=az+90 # assumed dip direction
#            #print('dd: ',dd[j]
#            if assdd>360:
#                assdd=assdd-360
#            #print('calc dd: ', assdd   ) 
#                if abs(assdd-dd[j])>=90:
#                    lon[j]=list(reversed(lon[j]))
#                    lat[j]=list(reversed(lat[j]))   
#                    #print('reversed'
#                
#        for j in range(len(lat)):
#            for k in range(len(lat[j])):
#                if k > (len(lat[j])-2):
#                    continue
#                lat1=lat[j][k]
#                lon1=lon[j][k]
#                lat2=lat[j][k+1]
#                lon2=lon[j][k+1]
#                az,baz,dist=g.inv(lon1,lat1,lon2,lat2,radians=False)
#                midlon,midlat,baz=g.fwd(lon1,lat1,az,dist/2,radians=False)
#                # SET FAULT PARAMETERS
#                if (az>=0):
#                    strike=az
#                else:
#                    strike=360+az
#    
#            
#                # build subfaults parameters using strikefactor and dipfactor
#                nstrike=int(round((dist/1000)/strikefactor))
#           	if (nstrike % 2 == 0):
#                    nstrike=nstrike+1
#                dx_strike=(dist/1000)/nstrike # in km
#                num_updip=0 # always zero because hypo is near surface
#                num_downdip=int(round(sd[j]/dipfactor))
#                dx_dip=(sd[j])/num_downdip
#                hypo=np.array([midlon,midlat,dx_dip/2])            
#                risetime=math.sqrt(dx_strike*dx_strike+dx_dip*dx_dip)/VR
#            
#                fout=home+'seg'+repr(j)+repr(k)+'.file' # name of output fault file
#                forward.makefault(fout,strike,dip[j],nstrike,dx_dip,dx_strike,hypo,num_updip,num_downdip,risetime)
#                
#        filelist=[]
#        filelist += [each for each in os.listdir(home) if each.endswith('.file')]
#        
#        outfile=open(home+directories[i]+'.fault','w')
#        outfile2=open(home+directories[i]+'.xy','w')
#        lastind=0
#        for file in filelist:
#            f=np.genfromtxt(home+file)   
#            for k in range(len(f)):
#                out='%i\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\n' % (f[k,0]+lastind,f[k,1],f[k,2],f[k,3],f[k,4],f[k,5],f[k,6],f[k,7],f[k,8],f[k,9])
#                outfile.write(out)
#                outfile2.write('%.6f\t%.6f\n' % (f[k,1],f[k,2]))
#                if (k==(len(f)-1)):
#                    lastind=f[-1,0]+lastind
#        print(lastind,' subfaults')
#        outfile.close()
#        outfile2.close()
        
def bssa14(M,Rjb,Vs30,U=0,RS=0,NS=0,SS=1,Z1=None,intensity_measure='SA'):
    '''
    Calculate ground motion intensity using the BSSA14 GMPE
    Note: This function will give the entire SA spectrum at a range of periods.
    
    Parameters:
        M - Moment magnitude
        Rjb - Distance to surface projection of fault in km
        U - is 1 if unspecified faulting style
        RS - is 1 if reverse faulting
        NS - is 1 if normal fault
        Vs30 - Vs30 in m/s
        Z1 - Depth to Vs=_1km/s, if unknown use Z1=None
        
    Returns:
        Y - the desired ground motion intensity, PGA in g, SA in g, or PGV in cm/s
        sigma - standard deviation
        periods - for SA, function also returns periods as first output
        
    Notes: 
        For strike slip faulting (default) set U=NS=RS=0
    '''

    from numpy import log,exp
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=float(Vs30)
    Rjb=float(Rjb)
    M=float(M)
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA':
        periods=[0.01,0.02,0.03,0.05,0.08,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.75,1.00,1.50,2.00,3.00,4.00,5.00,7.50,10.00]
        coefficients=[[0.4534,0.4916,0.2519,0.4599,1.421,0.04932,-0.1659,5.5,-1.13400,0.19160,-0.00809,4.5,1,4.5,0.000000,0.002820,-0.002440,-0.6037,1500.20,760,0,0.1,-0.1483,-0.00701,-9.9,-9.9,111.670,270.000,0.096,0.070,225,300,0.6980,0.4990,0.4020,0.3450],
                        [0.48598,0.52359,0.29707,0.48875,1.4331,0.053388,-0.16561,5.5,-1.13940,0.18962,-0.00807,4.5,1,4.5,0.000000,0.002780,-0.002340,-0.5739,1500.36,760,0,0.1,-0.1471,-0.00728,-9.9,-9.9,113.100,270.000,0.092,0.030,225,300,0.7020,0.5020,0.4090,0.3460],
                        [0.56916,0.6092,0.40391,0.55783,1.4261,0.061444,-0.1669,5.5,-1.14210,0.18842,-0.00834,4.5,1,4.49,0.000000,0.002760,-0.002170,-0.5341,1502.95,760,0,0.1,-0.1549,-0.00735,-9.9,-9.9,112.130,270.000,0.081,0.029,225,300,0.7210,0.5140,0.4450,0.3640],
                        [0.75436,0.79905,0.60652,0.72726,1.3974,0.067357,-0.18082,5.5,-1.11590,0.18709,-0.00982,4.5,1,4.2,0.000000,0.002960,-0.001990,-0.4580,1501.42,760,0,0.1,-0.1920,-0.00647,-9.9,-9.9,97.930,270.000,0.063,0.030,225,300,0.7530,0.5320,0.5030,0.4260],
                        [0.96447,1.0077,0.77678,0.9563,1.4174,0.073549,-0.19665,5.5,-1.08310,0.18225,-0.01058,4.5,1,4.04,0.000000,0.002960,-0.002160,-0.4441,1494.00,760,0,0.1,-0.2350,-0.00573,-9.9,-9.9,85.990,270.040,0.064,0.022,225,300,0.7450,0.5420,0.4740,0.4660],
                        [1.1268,1.1669,0.8871,1.1454,1.4293,0.055231,-0.19838,5.54,-1.06520,0.17203,-0.01020,4.5,1,4.13,0.000000,0.002880,-0.002440,-0.4872,1479.12,760,0,0.1,-0.2492,-0.00560,-9.9,-9.9,79.590,270.090,0.087,0.014,225,300,0.7280,0.5410,0.4150,0.4580],
                        [1.3095,1.3481,1.0648,1.3324,1.2844,-0.042065,-0.18234,5.74,-1.05320,0.15401,-0.00898,4.5,1,4.39,0.000000,0.002790,-0.002710,-0.5796,1442.85,760,0,0.1,-0.2571,-0.00585,-9.9,-9.9,81.330,270.160,0.120,0.015,225,300,0.7200,0.5370,0.3540,0.3880],
                        [1.3255,1.359,1.122,1.3414,1.1349,-0.11096,-0.15852,5.92,-1.06070,0.14489,-0.00772,4.5,1,4.61,0.000000,0.002610,-0.002970,-0.6876,1392.61,760,0,0.1,-0.2466,-0.00614,-9.9,-9.9,90.910,270.000,0.136,0.045,225,300,0.7110,0.5390,0.3440,0.3090],
                        [1.2766,1.3017,1.0828,1.3052,1.0166,-0.16213,-0.12784,6.05,-1.07730,0.13925,-0.00652,4.5,1,4.78,0.000000,0.002440,-0.003140,-0.7718,1356.21,760,0,0.1,-0.2357,-0.00644,-9.9,-9.9,97.040,269.450,0.141,0.055,225,300,0.6980,0.5470,0.3500,0.2660],
                        [1.2217,1.2401,1.0246,1.2653,0.95676,-0.1959,-0.092855,6.14,-1.09480,0.13388,-0.00548,4.5,1,4.93,0.000000,0.002200,-0.003300,-0.8417,1308.47,760,0,0.1,-0.2191,-0.00670,-9.9,-9.9,103.150,268.590,0.138,0.050,225,300,0.6750,0.5610,0.3630,0.2290],
                        [1.1046,1.1214,0.89765,1.1552,0.96766,-0.22608,-0.023189,6.2,-1.12430,0.12512,-0.00405,4.5,1,5.16,0.000000,0.002110,-0.003210,-0.9109,1252.66,760,0,0.1,-0.1958,-0.00713,-9.9,-9.9,106.020,266.540,0.122,0.049,225,300,0.6430,0.5800,0.3810,0.2100],
                        [0.96991,0.99106,0.7615,1.012,1.0384,-0.23522,0.029119,6.2,-1.14590,0.12015,-0.00322,4.5,1,5.34,0.000000,0.002350,-0.002910,-0.9693,1203.91,760,0,0.1,-0.1750,-0.00744,-9.9,-9.9,105.540,265.000,0.109,0.060,225,300,0.6150,0.5990,0.4100,0.2240],
                        [0.66903,0.69737,0.47523,0.69173,1.2871,-0.21591,0.10829,6.2,-1.17770,0.11054,-0.00193,4.5,1,5.6,0.000000,0.002690,-0.002530,-1.0154,1147.59,760,0,0.1,-0.1387,-0.00812,0.092,0.059,108.390,266.510,0.100,0.070,225,300,0.5810,0.6220,0.4570,0.2660],
                        [0.3932,0.4218,0.207,0.4124,1.5004,-0.18983,0.17895,6.2,-1.19300,0.10248,-0.00121,4.5,1,5.74,0.000000,0.002920,-0.002090,-1.0500,1109.95,760,0,0.1,-0.1052,-0.00844,0.367,0.208,116.390,270.000,0.098,0.020,225,300,0.5530,0.6250,0.4980,0.2980],
                        [-0.14954,-0.11866,-0.3138,-0.1437,1.7622,-0.1467,0.33896,6.2,-1.20630,0.09645,-0.00037,4.5,1,6.18,0.000000,0.003040,-0.001520,-1.0454,1072.39,760,0,0.1,-0.0620,-0.00771,0.638,0.309,125.380,262.410,0.104,0.010,225,300,0.5320,0.6190,0.5250,0.3150],
                        [-0.58669,-0.55003,-0.71466,-0.60658,1.9152,-0.11237,0.44788,6.2,-1.21590,0.09636,0.00000,4.5,1,6.54,0.000000,0.002920,-0.001170,-1.0392,1009.49,760,0,0.1,-0.0361,-0.00479,0.871,0.382,130.370,240.140,0.105,0.008,225,300,0.5260,0.6180,0.5320,0.3290],
                        [-1.1898,-1.142,-1.23,-1.2664,2.1323,-0.04332,0.62694,6.2,-1.21790,0.09764,0.00000,4.5,1,6.93,0.000000,0.002620,-0.001190,-1.0112,922.43,760,0,0.1,-0.0136,-0.00183,1.135,0.516,130.360,195.000,0.088,0.000,225,300,0.5340,0.6190,0.5370,0.3440],
                        [-1.6388,-1.5748,-1.6673,-1.7516,2.204,-0.014642,0.76303,6.2,-1.21620,0.10218,-0.00005,4.5,1,7.32,0.000000,0.002610,-0.001080,-0.9694,844.48,760,0,0.1,-0.0032,-0.00152,1.271,0.629,129.490,199.450,0.070,0.000,225,300,0.5360,0.6160,0.5430,0.3490],
                        [-1.966,-1.8882,-2.0245,-2.0928,2.2299,-0.014855,0.87314,6.2,-1.21890,0.10353,0.00000,4.5,1,7.78,0.000000,0.002600,-0.000570,-0.9195,793.13,760,0,0.1,-0.0003,-0.00144,1.329,0.738,130.220,230.000,0.061,0.000,225,300,0.5280,0.6220,0.5320,0.3350],
                        [-2.5865,-2.4874,-2.8176,-2.6854,2.1187,-0.081606,1.0121,6.2,-1.25430,0.12507,0.00000,4.5,1,9.48,0.000000,0.002600,0.000380,-0.7766,771.01,760,0,0.1,-0.0001,-0.00137,1.329,0.809,130.720,250.390,0.058,0.000,225,300,0.5120,0.6340,0.5110,0.2700],
                        [-3.0702,-2.9537,-3.3776,-3.1726,1.8837,-0.15096,1.0651,6.2,-1.32530,0.15183,0.00000,4.5,1,9.66,0.000000,0.003030,0.001490,-0.6558,775.00,760,0,0.1,0.0000,-0.00136,1.183,0.703,130.000,210.000,0.060,0.000,225,300,0.5100,0.6040,0.4870,0.2390]]   
    else:
        print('ERROR: Unknown intensity measure')
        #return
        
    if intensity_measure.upper()=='SA':
        YSA=[]
        sigmaSA=[]
        for i in range(len(periods)):
            #Assign each coefficient
            e0 = coefficients[i][0]
            e1 = coefficients[i][1]
            e2 = coefficients[i][2]
            e3 = coefficients[i][3]
            e4 = coefficients[i][4]
            e5 = coefficients[i][5]
            e6 = coefficients[i][6]
            Mh = coefficients[i][7]
            c1 = coefficients[i][8]
            c2 = coefficients[i][9]
            c3 = coefficients[i][10]
            Mref = coefficients[i][11]
            Rref = coefficients[i][12]
            h = coefficients[i][13]
            Dc3 = coefficients[i][14]
            Dc3chtur = coefficients[i][15]
            Dc3jpit = coefficients[i][16]
            C = coefficients[i][17]
            Vc = coefficients[i][18]
            Vref = coefficients[i][19]
            f1 = coefficients[i][20]
            f3 = coefficients[i][21]
            f4 = coefficients[i][22]
            f5 = coefficients[i][23]
            f6 = coefficients[i][24]
            f7 = coefficients[i][25]
            R1 = coefficients[i][26]
            R2 = coefficients[i][27]
            Dfr = coefficients[i][28]
            Dfv = coefficients[i][29]
            V1 = coefficients[i][30]
            V2 = coefficients[i][31]
            phi1 = coefficients[i][32]
            phi2 = coefficients[i][33]
            tau1 = coefficients[i][34]
            tau2 = coefficients[i][35]
        
            # Magnitude Scaling Term
            if NS == 0 and RS == 0 and U == 0:
                SS = 1
            else:
                SS = 0
        
            # Hinge magnitude term
            if M <= Mh:
                fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
            else:
                fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)  
            
            #Disance term
            R = (Rjb**2 + h**2)**0.5
        
            # Region term
            fp = (c1 + c2*(M - Mref))*log(R/Rref) + (c3 + Dc3)*(R - Rref)
        
            #Linear Site Term
            if Vs30 <= Vc:
                flin = C*log(Vs30/Vref)
            else:
                flin = C*log(Vc / Vref)
        
            #Nonlinear Site Term
            if Vs30 < 760:
                minV = Vs30
            else:
                minV = 760
        
            #Combine terms
            PGAr=bssa14_calc_one_station(M, Rjb, U, RS, NS,intensity_measure=intensity_measure)
            
            f2 = f4*((exp(f5*(minV - 360))) - exp(f5*(760 - 360)))
            fnl = f1 + f2*log((PGAr + f3)/f3)
            fnl = f1 + (f4*((exp(f5*(minV - 360)))-exp(f5*(760 - 360))))*log((PGAr + f3)/f3)
        
            #Basin Depth Term
            mz1 = exp(-7.15/4*log((Vs30**4 + 570.94**4)/(1360**4 + 570.94**4)))/1000
        
            #Final correction
            if Z1 == None:
                dz1 = 0
            else:
                dz1 = Z1 - mz1
            fz1 = 0
        
            if Z1 == None:
                fz1 = 0
            else:
                fz1 = fz1
        
            #Site Term
            fs = flin + fnl #in ln units
        
            #Model Prediction in ln units
            Y = exp(fm + fp + fs + fz1)            
            #Standard deviation
            sigma=bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure=intensity_measure)
            # build array
            YSA.append(Y)
            sigmaSA.append(sigma)
            
        Y=YSA
        sigma=sigmaSA
        return periods,Y,sigma
        
    else:
        #Assign each coefficient
        e0 = coefficients[0]
        e1 = coefficients[1]
        e2 = coefficients[2]
        e3 = coefficients[3]
        e4 = coefficients[4]
        e5 = coefficients[5]
        e6 = coefficients[6]
        Mh = coefficients[7]
        c1 = coefficients[8]
        c2 = coefficients[9]
        c3 = coefficients[10]
        Mref = coefficients[11]
        Rref = coefficients[12]
        h = coefficients[13]
        Dc3 = coefficients[14]
        Dc3chtur = coefficients[15]
        Dc3jpit = coefficients[16]
        C = coefficients[17]
        Vc = coefficients[18]
        Vref = coefficients[19]
        f1 = coefficients[20]
        f3 = coefficients[21]
        f4 = coefficients[22]
        f5 = coefficients[23]
        f6 = coefficients[24]
        f7 = coefficients[25]
        R1 = coefficients[26]
        R2 = coefficients[27]
        Dfr = coefficients[28]
        Dfv = coefficients[29]
        V1 = coefficients[30]
        V2 = coefficients[31]
        phi1 = coefficients[32]
        phi2 = coefficients[33]
        tau1 = coefficients[34]
        tau2 = coefficients[35]
    
        # Magnitude Scaling Term
        if NS == 0 and RS == 0 and U == 0:
            SS = 1
        else:
            SS = 0
    
        # Hinge magnitude term
        if M <= Mh:
            fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
        else:
            fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)  
        
        #Disance term
        R = (Rjb**2 + h**2)**0.5
    
        # Region term
        fp = (c1 + c2*(M - Mref))*log(R/Rref) + (c3 + Dc3)*(R - Rref)
    
        #Linear Site Term
        if Vs30 <= Vc:
            flin = C*log(Vs30/Vref)
        else:
            flin = C*log(Vc / Vref)
    
        #Nonlinear Site Term
        if Vs30 < 760:
            minV = Vs30
        else:
            minV = 760
    
        #Combine terms
        PGAr=bssa14_calc_one_station(M, Rjb, U, RS, NS,intensity_measure=intensity_measure)
        
        f2 = f4*((exp(f5*(minV - 360))) - exp(f5*(760 - 360)))
        fnl = f1 + f2*log((PGAr + f3)/f3)
        fnl = f1 + (f4*((exp(f5*(minV - 360)))-exp(f5*(760 - 360))))*log((PGAr + f3)/f3)
    
        #Basin Depth Term
        mz1 = exp(-7.15/4*log((Vs30**4 + 570.94**4)/(1360**4 + 570.94**4)))/1000
    
        #Final correction
        if Z1 == None:
            dz1 = 0
        else:
            dz1 = Z1 - mz1
    
        fz1 = 0
    
        if Z1 == None:
            fz1 = 0
        else:
            fz1 = fz1
    
        #Site Term
        fs = flin + fnl #in ln units
    
        #Model Prediction in ln units        
        Y = exp(fm + fp + fs + fz1)
        
        #Standard deviation
        sigma=bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure=intensity_measure)
        
        return Y,sigma

def bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure='PGA'):
    '''
    Get GMPE standard deviation
    '''
    from numpy import log,exp,sqrt,array,ones,where,zeros
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=float(Vs30)
    Rjb=float(Rjb)
    M=float(M)
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA':
        coefficients=[-3.0702,-2.9537,-3.3776,-3.1726,1.8837,-0.15096,1.0651,6.2,-1.32530,
                        0.15183,0.00000,4.5,1,9.66,0.000000,0.003030,0.001490,-0.6558,775.00,
                        760,0,0.1,0.0000,-0.00136,1.183,0.703,130.000,210.000,0.060,0.000,225,
                        300,0.5100,0.6040,0.4870,0.2390]
    else:
        print('ERROR: Unknown intensity measure')
        #return

    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]


    if M <= 4.5:
        tauM = tau1
    elif M > 4.5 and M < 5.5:
        tauM = tau1 + (tau2 - tau1) * (M - 4.5)
    else:
        tauM = tau2
    
    if M <= 4.5:
        phiM = phi1
    elif M > 4.5 and M < 5.5:
        phiM = phi1 + (phi2 - phi1) * (M - 4.5)
    else:
        phiM = phi2

    
    if Rjb <= R1:
        phiMR = phiM
    elif Rjb > R1 and Rjb <= R2:
        phiMR = phiM + Dfr * (log(Rjb / R1) / (log(R2 / R1)))
    else:
        phiMR = phiM + Dfr
    
    if Vs30 >= V2:
        phi = phiMR
    elif Vs30 >= V1 and Vs30 <= V2:
        phi = phiMR - Dfv * (log(V2 / Vs30) / (log(V2 / V1)))
    else:
        phi = phiMR - Dfv
    
    #Model Prediction in ln units
    sigma = (tauM**2 + phi**2)**0.5
    
    return sigma

def bssa14_calc_one_station(M, Rjb, U, RS, NS,intensity_measure):
    '''
    Calculate reference PGA
    '''
    
    from numpy import log,exp,ones
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV, T(s)=10 s for SA
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA':
        coefficients=[-3.0702,-2.9537,-3.3776,-3.1726,1.8837,-0.15096,1.0651,6.2,-1.32530,
                        0.15183,0.00000,4.5,1,9.66,0.000000,0.003030,0.001490,-0.6558,775.00,
                        760,0,0.1,0.0000,-0.00136,1.183,0.703,130.000,210.000,0.060,0.000,225,
                        300,0.5100,0.6040,0.4870,0.2390]
    else:
        print('ERROR: Unknown intensity measure')
        #return
        
    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    
    if NS == 0 and RS == 0 and U == 0:
        SS = 1
    else:
        SS = 0

    if M <= Mh:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
    else:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)
    
    R = (Rjb**2 + h**2)**0.5
    
    #region term
    fp = (c1 + c2 * (M - Mref)) * log(R / Rref) + (c3 + Dc3) * (R - Rref)
    
    #Calculate PGAr
    PGAr = exp(fm + fp)
    
    return PGAr
    
def bssa14_scalar(M,Rjb,Vs30,U=0,RS=0,NS=0,SS=1,Z1=None,intensity_measure='SA'):
    '''
    Calculate ground motion intensity using the BSSA14 GMPE
    
    Parameters:
        M - Moment magnitude
        Rjb - Distance to surface projection of fault in km
        U - is 1 if unspecified faulting style
        RS - is 1 if reverse faulting
        NS - is 1 if normal fault
        Vs30 - Vs30 in m/s
        Z1 - Depth to Vs=_1km/s, if unknown use Z1=None
        
    Returns:
        Y - the desired ground motion intensity, PGA in g, SA in g, or PGV in cm/s
        
    Notes: 
        For strike slip faulting (default) set U=NS=RS=0
    '''

    from numpy import log,exp
    
    # GMPE coefficients from the PEER spreadsheet: 
    # http://peer.berkeley.edu/ngawest2/wp-content/uploads/2016/02/NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.zip
    # in the "BSSA14_Coeffs sheet, T(s)=0 corresponds to PGA, T(s)=-1 is PGV
    
    #Convert input to floats
    Vs30=float(Vs30)
    Rjb=float(Rjb)
    M=float(M)
    
    if intensity_measure.upper()=='PGA':
        coefficients=[0.4473,0.4856,0.2459,0.4539,1.431,0.05053,-0.1662,5.5,-1.13400,0.19170,-0.00809,4.5,
                        1.,4.5,0.000000,0.002860,-0.002550,-0.6000,1500.00,760,0.,0.1,-0.1500,-0.00701,-9.900,
                        -9.900,110.000,270.000,0.100,0.070,225.,300.,0.6950,0.4950,0.3980,0.3480]
    elif intensity_measure.upper()=='PGV':
        coefficients=[5.037,5.078,4.849,5.033,1.073,-0.1536,0.2252,6.2,-1.24300,0.14890,-0.00344,4.5,1.,5.3,
                        0.000000,0.004350,-0.000330,-0.8400,1300.00,760,0.,0.1,-0.1000,-0.00844,-9.900,-9.900,
                        105.000,272.000,0.082,0.080,225.,300.,0.6440,0.5520,0.4010,0.3460]
    elif intensity_measure.upper()=='SA':
        coefficients=[-3.0702,-2.9537,-3.3776,-3.1726,1.8837,-0.15096,1.0651,6.2,-1.32530,
                        0.15183,0.00000,4.5,1,9.66,0.000000,0.003030,0.001490,-0.6558,775.00,
                        760,0,0.1,0.0000,-0.00136,1.183,0.703,130.000,210.000,0.060,0.000,225,
                        300,0.5100,0.6040,0.4870,0.2390]
    else:
        print('ERROR: Unknown intensity measure')
        #return
        
    #Assign each coefficient
    e0 = coefficients[0]
    e1 = coefficients[1]
    e2 = coefficients[2]
    e3 = coefficients[3]
    e4 = coefficients[4]
    e5 = coefficients[5]
    e6 = coefficients[6]
    Mh = coefficients[7]
    c1 = coefficients[8]
    c2 = coefficients[9]
    c3 = coefficients[10]
    Mref = coefficients[11]
    Rref = coefficients[12]
    h = coefficients[13]
    Dc3 = coefficients[14]
    Dc3chtur = coefficients[15]
    Dc3jpit = coefficients[16]
    C = coefficients[17]
    Vc = coefficients[18]
    Vref = coefficients[19]
    f1 = coefficients[20]
    f3 = coefficients[21]
    f4 = coefficients[22]
    f5 = coefficients[23]
    f6 = coefficients[24]
    f7 = coefficients[25]
    R1 = coefficients[26]
    R2 = coefficients[27]
    Dfr = coefficients[28]
    Dfv = coefficients[29]
    V1 = coefficients[30]
    V2 = coefficients[31]
    phi1 = coefficients[32]
    phi2 = coefficients[33]
    tau1 = coefficients[34]
    tau2 = coefficients[35]

    # Magnitude Scaling Term
    if NS == 0 and RS == 0 and U == 0:
        SS = 1
    else:
        SS = 0

    # Hinge magnitude term
    if M <= Mh:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e4*(M - Mh) + e5*(M - Mh)**2
    else:
        fm = e0*U + e1*SS + e2*NS + e3*RS + e6*(M - Mh)  
    
    #Disance term
    R = (Rjb**2 + h**2)**0.5

    # Region term
    fp = (c1 + c2*(M - Mref))*log(R/Rref) + (c3 + Dc3)*(R - Rref)

    #Linear Site Term
    if Vs30 <= Vc:
        flin = C*log(Vs30/Vref)
    else:
        flin = C*log(Vc / Vref)

    #Nonlinear Site Term
    if Vs30 < 760:
        minV = Vs30
    else:
        minV = 760

    #Combine terms
    PGAr=bssa14_calc_one_station(M, Rjb, U, RS, NS,intensity_measure=intensity_measure)
    
    f2 = f4*((exp(f5*(minV - 360))) - exp(f5*(760 - 360)))
    fnl = f1 + f2*log((PGAr + f3)/f3)
    fnl = f1 + (f4*((exp(f5*(minV - 360)))-exp(f5*(760 - 360))))*log((PGAr + f3)/f3)

    #Basin Depth Term
    mz1 = exp(-7.15/4*log((Vs30**4 + 570.94**4)/(1360**4 + 570.94**4)))/1000

    #Final correction
    if Z1 == None:
        dz1 = 0
    else:
        dz1 = Z1 - mz1

    fz1 = 0

    if Z1 == None:
        fz1 = 0
    else:
        fz1 = fz1

    #Site Term
    fs = flin + fnl #in ln units

    #Model Prediction in ln units
    Y = exp(fm + fp + fs + fz1)
    
    #Standard deviation
    sigma=bssa14_stdev_one_station(M,Rjb,Vs30,intensity_measure=intensity_measure)
    
    return Y,sigma
        
def calc_SA(home,project_name,run_name,GF_list,fault_name,rupture_list='ruptures.list',slip_type='SS',Vs30=760):
    '''
    Calculate response spectra from waveforms as well as the spectral acceleration predicted by bssa14 GMPE.
    Write to files in waveforms directories.
    
    Output used for validation.
    
    Christine J. Ruhl, August 2016
    '''
    import numpy as np
    from numpy import genfromtxt
    from obspy import read
    
    # start of function
    frequency_vector=np.array([0.1]) # only do 10 s period
    
    # Read summary file
    if project_name=='Cascadia':
        sta=genfromtxt(home+project_name+'/'+GF_list,usecols=0,dtype='S')
        slon=genfromtxt(home+project_name+'/'+GF_list,usecols=1,dtype='float')
        slat=genfromtxt(home+project_name+'/'+GF_list,usecols=2,dtype='float')
        rupt=genfromtxt(home+project_name+'/'+rupture_list,usecols=0,dtype='S')
    else:
        sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
        slon=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=1,dtype='float')
        slat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=2,dtype='float')
        rupt=genfromtxt(home+project_name+'/data/'+rupture_list,usecols=0,dtype='S')
    for j in range(len(rupt)):
        print('calculating SA for rupture '+rupt[j].split('.')[1] )
        # open summary file for writing
        # save in home+project_name+'/output/waveforms/_analysis.'+run_name+'.'+rupt[j].split('.')[1]+'.txt'
        if project_name=='Cascadia':
            sf=open(home+project_name+'/'+run_name+'.'+rupt[j].split('.')[1]+'/_analysis.'+run_name+'.'+rupt[j].split('.')[1]+'.txt','w')
        else:
            sf=open(home+project_name+'/output/waveforms/'+run_name+'.'+rupt[j].split('.')[1]+'/_analysis.'+run_name+'.'+rupt[j].split('.')[1]+'.txt','w')
        line_out='# Station, lon, lat, SAobsEAST [g], SAobsNORTH [g], SAobsRotD50 [g], SAcalc [g], sigma [g], Rjb [km]\n'
            
        for k in range(len(sta)):            
            ### Calculate new response spectra using pyrotd code responseSpectrum
            ### responseSpectrum(timeStep, accelTs, oscFreqs, oscDamping=0.05)
            if project_name=='Cascadia':
                e=read(home+project_name+'/'+run_name+'.'+rupt[j].split('.')[1]+'/'+sta[k]+'.LYE.sac')
                n=read(home+project_name+'/'+run_name+'.'+rupt[j].split('.')[1]+'/'+sta[k]+'.LYN.sac')
            else:  
                e=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt[j].split('.')[1]+'/'+sta[k]+'.LYE.sac')
                n=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt[j].split('.')[1]+'/'+sta[k]+'.LYN.sac')
                
            # diff twice to get acceleration
            e.differentiate()
            e.differentiate()
            # calculate response
            respe=responseSpectrum(e[0].stats.sampling_rate,e[0].data/9.81,frequency_vector)          
            # diff twice to get acceleration
            n.differentiate()
            n.differentiate()
            # calculate response
            respn=responseSpectrum(n[0].stats.sampling_rate,n[0].data/9.81,frequency_vector)
            
            # average them                
            resp=rotatedResponseSpectrum(n[0].stats.sampling_rate, n[0].data/9.81, e[0].data/9.81, frequency_vector)
  
            # get the magnitude of each rupture from the log file
            if project_name=='Cascadia':
                logfile=home+project_name+'/'+run_name+'.'+rupt[j].split('.')[1]+'/_'+run_name+'.'+rupt[j].split('.')[1]+'.log'
            else:
                logfile=home+project_name+'/output/ruptures/'+run_name+'.'+rupt[j].split('.')[1]+'.log'
            f=open(logfile,'r')
            loop_go=True
            while loop_go:
                line=f.readline()  
                if 'Actual magnitude' in line:
                    Mw=float(line.split(':')[-1].split(' ')[-1])
                    break
                    
            # calculate Rjb
            if project_name=='Cascadia':
                ruptfile=home+project_name+'/'+run_name+'.'+rupt[j].split('.')[1]+'/_'+run_name+'.'+rupt[j].split('.')[1]+'.rupt'
                #print(ruptfile)
            else:
                ruptfile=home+project_name+'/output/ruptures/'+run_name+'.'+rupt[j].split('.')[1]+'.rupt'
                
            Rjb=calc_rjb(slon[k],slat[k],ruptfile)
            
            # set Vs30
            Vs30=np.array(Vs30) # keep constant?
            # calculate SA based on slip_type
            if slip_type=='SS':
                sa,std=bssa14_scalar(Mw,Rjb,Vs30,U=0,RS=0,NS=0,SS=1,Z1=None,intensity_measure='SA')
            elif slip_type=='NS':
                sa,std=bssa14_scalar(Mw,Rjb,Vs30,U=0,RS=0,NS=1,SS=0,Z1=None,intensity_measure='SA')
            elif slip_type=='RS':
                sa,std=bssa14_scalar(Mw,Rjb,Vs30,U=0,RS=1,NS=0,SS=0,Z1=None,intensity_measure='SA')       
            
            # print(to summary file in rupture folder)
            # convert response spectral acceleration to g
            line_out+='%s\t%10.4f\t%10.4f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n' % (sta[k],slon[k],slat[k],respe[0],respn[0],resp[0][0][0],sa,std,Rjb)

        sf.write(line_out)
        sf.close()

def one_event_sa_misfit(home,project_name,run_name,rupt):
    '''
    Plot SA misfit from one event.
    
    Christine J. Ruhl, August 2016
    '''
    
    import numpy as np
    from numpy import genfromtxt
    from matplotlib import pyplot as plt
    
    # Read summary file
    analysis_file=home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/_analysis.'+run_name+'.'+rupt+'.txt'
    SArotD50=genfromtxt(analysis_file,usecols=[5])
    SAcalc=genfromtxt(analysis_file,usecols=[6])
    sigma=genfromtxt(analysis_file,usecols=[7])
    rjb=genfromtxt(analysis_file,usecols=[8])
    
    # get the magnitude of each rupture from the log file
    logfile=home+project_name+'/output/ruptures/'+run_name+'.'+rupt+'.log'
    f=open(logfile,'r')
    loop_go=True
    while loop_go:
        line=f.readline()  
        if 'Actual magnitude' in line:
            Mw=float(line.split(':')[-1].split(' ')[-1])
            break

    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,sharey=False)
    ax1.loglog(rjb,SArotD50,'k.')
    ax1.set_ylabel('SA (g)'), ax1.set_xlabel('Rjb (km)')
    
    # sort by Rjb
    ind=[i[0] for i in sorted(enumerate(rjb), key=lambda x:x[1])]
    ax1.loglog(rjb[ind],SAcalc[ind],'k-')
    ax1.loglog(rjb[ind],(np.exp(np.log(SAcalc[ind])-np.log(sigma[ind]))),'k--')
    ax1.loglog(rjb[ind],(np.exp(np.log(SAcalc[ind])+np.log(sigma[ind]))),'k--')
    ax1.set_title(run_name+' '+rupt+' Mw '+str(Mw))
    
    # plot residuals
    ax2.loglog(rjb[ind],SArotD50[ind]/SAcalc[ind],'k.')
    ax2.set_ylabel('Obs/Calc Misfit'), ax2.set_xlabel('Rjb (km)')
    
    plt.show() 

def SA_2D_misfit_and_GOF(home,project_name,rupture_list='ruptures.list',Mw_lims=[7.8,9.5],dist_lims=[10,1000],misfit_lims=[-3,3],GOF_lims=[0,2],n_mag_bins=10,n_dist_bins=10):
    '''
    Plot spectral acceleration misfit as a function of both distance and magnitude
    '''
    
    from numpy import genfromtxt,array,zeros,logspace,linspace,r_,log,genfromtxt,log10,ones,where,arange,mean
    from matplotlib import pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import MaxNLocator
    
    ruptures=genfromtxt(home+project_name+'/data/'+rupture_list,dtype='S')
    SAobs_all=array([])
    rjb_all=array([])
    Mw_all=array([])
    SAcalc_all=array([])
    
    for k in range(len(ruptures)):
        run_name=ruptures[k].split('.')[0]
        run_number=ruptures[k].split('.')[1]
        # Read analysis file
        analysis_file=home+project_name+'/output/waveforms/'+run_name+'.'+run_number+'/_analysis.'+run_name+'.'+run_number+'.txt'
        SAobs=genfromtxt(analysis_file,usecols=[5])      
        SAcalc=genfromtxt(analysis_file,usecols=[6])
        rjb=genfromtxt(analysis_file,usecols=[8])
        
        # get the magnitude of each rupture from the log file
        logfile=home+project_name+'/output/ruptures/'+run_name+'.'+run_number+'.log'
        f=open(logfile,'r')
        loop_go=True
        while loop_go:
            line=f.readline()  
            if 'Actual magnitude' in line:
                Mw=float(line.split(':')[-1].split(' ')[-1])
                break
           
        #Concatente to output variables
        SAobs_all=r_[SAobs_all,SAobs]
        rjb_all=r_[rjb_all,rjb]
        Mw_all=r_[Mw_all,Mw*ones(len(rjb))]                   
        SAcalc_all=r_[SAcalc_all,SAcalc] 
                    
    #Get misfits             
    misfit=log(SAobs_all/SAcalc_all)                    

    #remove extrneous magnitudes
    i=where((Mw_all>=Mw_lims[0]) & (Mw_all<=Mw_lims[1]))[0]
    Mw_all=Mw_all[i]
    misfit=misfit[i]
    rjb_all=rjb_all[i]
    
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
            i=where((Mw_all>=bin_edges_x[kx]) & (Mw_all<bin_edges_x[kx+1]) & (rjb_all>=bin_edges_y[ky]) & (rjb_all<bin_edges_y[ky+1]))[0]
            GOF_bin[kx,ky]=0.5*(abs(mean(misfit[i])))+0.5*mean(abs(misfit[i]))
            misfit_bin[kx,ky]=misfit[i].mean()
            frequency_bin[kx,ky]=len(misfit[i])
            bin_centers_x[kx]=bin_edges_x[kx+1]-((bin_edges_x[kx+1]-bin_edges_x[kx])/2)        
            bin_centers_y[ky]=bin_edges_y[ky+1]-((bin_edges_y[ky+1]-bin_edges_y[ky])/2) 
            
    plt.figure(figsize=(15,5.5))    
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
    plt.ylabel('Rjb Distance (km)')
    levs=ax.get_xticks()
    ax.set_xticklabels(levs,rotation=-55)
    #plt.annotate('(a)',xy=(7.9,680),fontsize=16)
    bbox_props = dict(boxstyle="round", fc="w")
    ax.text(7.9, 350, "(a)", ha="center", va="center", size=18,bbox=bbox_props)
    cb = plt.colorbar(PM, orientation='horizontal',pad = 0.02) 
    cb.set_label('Ln(10s SA Obs/Calc)')
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
    plt.show()    
    
## Functions borrowed from pyrotd
def oscillatorTimeSeries(freq, fourierAmp, oscFreq, oscDamping):
    '''
    Function borrowed from pyrotd.
    
    Compute the time series response of an oscillator.

    Parameters
    ----------
    freq: numpy.array
        frequency of the Fourier acceleration spectrum [Hz]
    fourierAmp: numpy.array
        Fourier acceleration spectrum [g-sec]
    oscFreq: float
        frequency of the oscillator [Hz]
    oscDamping: float
        damping of the oscillator [decimal]

    Returns
    -------
    response: numpy.array
        time series response of the oscillator
    '''
    
    from numpy import power,fft
    
    # Single-degree of freedom transfer function
    h = (-power(oscFreq, 2.)
            / ((power(freq, 2.) - power(oscFreq, 2.))
                - 2.j * oscDamping * oscFreq * freq))

    # Adjust the maximum frequency considered. The maximum frequency is 5
    # times the oscillator frequency. This provides that at the oscillator
    # frequency there are at least tenth samples per wavelength.
    n = len(fourierAmp)
    m = max(n, int(2. * oscFreq / freq[1]))
    scale = float(m) / float(n)

    # Scale factor is applied to correct the amplitude of the motion for the
    # change in number of points
    return scale * fft.irfft(fourierAmp * h, 2 * (m-1))

def peakResponse(resp):
    '''
    Function borrowed from pyrotd.
    
    Compute the maximum absolute value of a response.

    Parameters
    ----------
    resp: numpy.array
        time series of a response

    Returns
    -------
    peakResponse: float
        peak response
    '''
    from numpy import max,abs
    
    return max(abs(resp))

def rotatedResponseSpectrum(timeStep, accelA, accelB, oscFreqs, oscDamping=0.05,
        percentiles=[50], angles=np.arange(0, 180, step=1)):
    '''
    Function borrowed from pyrotd.
    
    Compute the response spectrum for a time series.

    Parameters
    ----------
    timeStep: float
        time step of the time series [s]
    accelA: numpy.array
        acceleration time series of the first motion [g]
    accelB: numpy.array
        acceleration time series of the second motion that is perpendicular to the first motion [g]
    oscFreqs: numpy.array
        natural frequency of the oscillators [Hz]
    oscDamping: float
        damping of the oscillator [decimal]. Default of 0.05 (i.e., 5%)
    percentiles: numpy.array
        percentiles to return. Default of [0, 50, 100],
    angles: numpy.array
        angles to which to compute the rotated time series. Default of
        np.arange(0, 180, step=1) (i.e., 0, 1, 2, .., 179).

    Returns
    -------
    oscResps: list(numpy.array)
        computed psuedo-spectral acceleartion [g] at each of the percentiles
    '''

    from numpy import array

    assert len(accelA) == len(accelB), 'Time series not equal lengths!'

    # Compute the Fourier amplitude spectra
    fourierAmps = [np.fft.rfft(accelA), np.fft.rfft(accelB)]
    freq = np.linspace(0, 1./(2 * timeStep), num=fourierAmps[0].size)

    values = []
    for i, oscFreq in enumerate(oscFreqs):
        # Compute the oscillator responses
        oscResps = [oscillatorTimeSeries(freq, fa, oscFreq, oscDamping)
                for fa in fourierAmps]

        # Compute the rotated values of the oscillator response
        vals,orients = rotatedPercentiles(oscResps[0], oscResps[1], angles, percentiles)
        values.append(vals[0])

    # Reorganzie the arrays grouping by the percentile
#    oscResps = [np.array([v[i] for v in values],
#        dtype=[('value', '<f8'), ('orientation', '<f8')]) for i in range(len(percentiles))]

    return array(values)

def rotatedPercentiles(accelA, accelB, angles, percentiles=[50]):
    '''
    Function borrowed from pyrotd.
    
    Compute the response spectrum for a time series.

    Parameters
    ----------
    accelA: numpy.array
        first time series
    accelB: numpy.array
        second time series that is perpendicular to the first
    angles: numpy.array
        angles to which to compute the rotated time series
    percentiles: numpy.array
        percentiles to return

    Returns
    -------
    values: numpy.array
        rotated values and orientations corresponding to the percentiles
    '''
    assert all(0 <= p <= 100 for p in percentiles), 'Invalid percentiles.'

    # Compute the response for each of the specified angles and sort this array
    # based on the response
    rotated = np.array(
            [(a, peakResponse(rotateTimeSeries(accelA, accelB, a))) for a in angles],
            dtype=[('angle', '<f8'), ('value', '<f8')])
    rotated.sort(order='value')

    # Interpolate the percentile from the values
    values = np.interp(percentiles,
            np.linspace(0, 100, len(angles)), rotated['value'])

    # Can only return the orientations for the minimum and maximum value as the
    # orientation is not unique (i.e., two values correspond to the 50%
    # percentile).
    orientationMap = {
            0 : rotated['angle'][0],
            100 : rotated['angle'][-1],
            }
    orientations = [orientationMap.get(p, np.nan) for p in percentiles] 

#    out=np.array(zip(values, orientations), dtype=[('value', '<f8'), ('orientation', '<f8')])

    return values,orientations
            
def rotateTimeSeries(foo, bar, angle):
    '''
    Function borrowed from pyrotd.
    
    Compute the rotated time series.

    Parameters
    ----------
    foo: numpy.array
        first time series
    bar: numpy.array
        second time series that is perpendicular to the first

    Returns
    -------
    foobar: numpy.array
        time series rotated by the specified angle
    '''

    angleRad = np.radians(angle)
    # Rotate the time series using a vector rotation
    return foo * np.cos(angleRad) + bar * np.sin(angleRad)
    
def responseSpectrum(timeStep, accelTs, oscFreqs, oscDamping=0.05):
    '''
    Function borrowed from pyrotd.
    
    Compute the response spectrum for a time series.

    Parameters
    ----------
    timeStep: float
        time step of the time series [s]
    accelTs: numpy.array
        acceleration time series [g]
    oscFreqs: numpy.array
        natural frequency of the oscillators [Hz]
    oscDamping: float
        damping of the oscillator [decimal]. Default of 0.05 (i.e., 5%)

    Returns
    -------
    oscResp: numpy.array
        computed psuedo-spectral acceleartion [g]
    '''
    fourierAmp = np.fft.rfft(accelTs)
    freq = np.linspace(0, 1./(2 * timeStep), num=fourierAmp.size)

    psa = [peakResponse(oscillatorTimeSeries(freq, fourierAmp, of, oscDamping))
            for of in oscFreqs]

    return np.array(psa)

def sac2rt_textfile(home,project_name,run_name,time_epi,rupture_list):
    '''
    Function to write sac waveforms to textfile for realtime (rt) simulation in G-larmS.
    
    Christine J. Ruhl, September 2016    
    '''
    from obspy import read
    import numpy as np
    
    # load ruptures
    if project_name=='Cascadia':
        ruptures=np.genfromtxt(home+project_name+'/'+rupture_list,usecols=0,dtype='S')
    else:
        ruptures=np.genfromtxt(home+project_name+'/data/'+rupture_list,usecols=0,dtype='S')
    
    for j in range(len(ruptures)):
        # get rupture number
        rupt=ruptures[j].split('.')[1]
        print('writing realtime textfiles for rupture '+rupt)
        if project_name=='Cascadia':
            # open summary file for reading
            sumfile=home+project_name+'/'+run_name+'.'+rupt+'/_'+run_name+'.'+rupt+'.offsets'
            # open new offsetf file for writing
            outfile1=home+project_name+'/'+run_name+'.'+rupt+'/_horiz.'+run_name+'.'+rupt+'.txt'
            outfile2=home+project_name+'/'+run_name+'.'+rupt+'/_vert.'+run_name+'.'+rupt+'.txt'
        else:
            # open summary file for reading
            sumfile=home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/_summary.'+run_name+'.'+rupt+'.txt'
            # open new offsetf file for writing
            outfile1=home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/_horiz.'+run_name+'.'+rupt+'.txt'
            outfile2=home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/_vert.'+run_name+'.'+rupt+'.txt'
            
        file1=open(outfile1,'w')
        file2=open(outfile2,'w')
        
        sta=np.genfromtxt(sumfile,usecols=0,dtype='S')
        lonlat=np.genfromtxt(sumfile,usecols=[1,2])
        
        for i in range(len(sta)):
            # read waveforms
            if project_name=='Cascadia':
                e=read(home+project_name+'/'+run_name+'.'+rupt+'/'+sta[i]+'.LYE.sac')
                n=read(home+project_name+'/'+run_name+'.'+rupt+'/'+sta[i]+'.LYN.sac')
                z=read(home+project_name+'/'+run_name+'.'+rupt+'/'+sta[i]+'.LYZ.sac')
            else:
                e=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[i]+'.LYE.sac')
                n=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[i]+'.LYN.sac')
                z=read(home+project_name+'/output/waveforms/'+run_name+'.'+rupt+'/'+sta[i]+'.LYZ.sac')
            delta=e[0].stats.starttime-time_epi
            
            # print(adjusted times)
            for j in range(len(e[0].data)):
                file1.write('%.6f %.6f %.6f %.6f %.6f %d %d %d %s\n' % (e[0].times()[j]+delta,
                    lonlat[i][0],lonlat[i][1],e[0].data[j],n[0].data[j],0,0,0,sta[i]))
                file2.write('%.6f %.6f %.6f %.6f %.6f %d %d %d %s\n' % (z[0].times()[j]+delta,
                    lonlat[i][0],lonlat[i][1],0,z[0].data[j],0,0,0,sta[i]))
                    
        file1.close()
        file2.close()
        sort_rt_textfile(outfile1)
        sort_rt_textfile(outfile2)

def sort_rt_textfile(outfile):
    '''
    Function to resort realtime textfile:
    
    Sort realtime textfile (called outfile) from:
        TIME1 lon lat e n esig nsig en-corr STA1
        TIME​2​ lon lat e n esig nsig en-corr STA1
        TIME3 lon lat e n esig nsig en-corr STA1
        TIME1 lon lat e n esig nsig en-corr STA2
        TIME2 lon lat e n esig nsig en-corr STA2
        TIME​3​ lon lat e n esig nsig en-corr STA​2​
    to:
        TIME1 lon lat e n esig nsig en-corr STA1
        TIME1 lon lat e n esig nsig en-corr STA2
        TIME​2​ lon lat e n esig nsig en-corr STA1
        TIME2 lon lat e n esig nsig en-corr STA2
        TIME3 lon lat e n esig nsig en-corr STA1
        TIME​3​ lon lat e n esig nsig en-corr STA​2​
        
    Christine J. Ruhl, September 2016
    '''
    # example outfile1 = '/Users/cruhl/Projects/fakequakes/Hayward_1Hz_newvel/output/waveforms/hayward.000149/_vert.hayward.000149.txt'
    lines=sorted(open(outfile).readlines(), key=lambda line: float(line.split(' ')[0]))
    f=open(outfile,'w')
    for i in range(len(lines)):
        f.write(lines[i])
    f.close()
