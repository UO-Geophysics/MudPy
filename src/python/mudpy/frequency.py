'''
D.Melgar 05/2014

Frequency domain analysis tools
'''

def allpsd(home,project_name,run_name,run_number,GF_list,d_or_s,v_or_d,decimate,lowpass):
    '''
    Compute PSDs for either all the synthetics fo a aprticular run or all the data
    '''
    
    from numpy import genfromtxt,where,log10,savez
    from obspy import read
    from mudpy.forward import lowpass as lfilter
    from mudpy.green import stdecimate 
    import nitime.algorithms as tsa
    
    #Decide what I'm going to work on
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    outpath=home+project_name+'/analysis/frequency/'
    if v_or_d.lower()=='d':
        kgf=0 #disp
        datasuffix='kdisp'
        synthsuffix='disp'
    elif v_or_d.lower()=='v':
        kgf=1 #disp
        datasuffix='kvel'
        synthsuffix='vel'
    if d_or_s.lower()=='d': #We're working on observed data
        path=datapath
        suffix=datasuffix
    else: #We're looking at syntehtics from a certain run
        path=synthpath
        suffix=synthsuffix
    i=where(gf[:,kgf]==1)[0]
    for k in range(len(i)):
        print('Working on '+sta[i[k]])
        if d_or_s.lower()=='d': #Read data
            n=read(path+sta[i[k]]+'.'+suffix+'.n')
            e=read(path+sta[i[k]]+'.'+suffix+'.e')
            u=read(path+sta[i[k]]+'.'+suffix+'.u')
            outname=sta[i[k]]+'.'+suffix+'.psd'
            if lowpass!=None:
                fsample=1./e[0].stats.delta
                e[0].data=lfilter(e[0].data,lowpass,fsample,10)
                n[0].data=lfilter(n[0].data,lowpass,fsample,10)
                u[0].data=lfilter(u[0].data,lowpass,fsample,10)
            if decimate!=None:
                n[0]=stdecimate(n[0],decimate)
                e[0]=stdecimate(e[0],decimate)
                u[0]=stdecimate(u[0],decimate)
        else: #Read synthetics
            n=read(path+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.n.sac')
            e=read(path+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.e.sac')
            u=read(path+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.u.sac')
            outname=run_name+'.'+run_number+'.'+sta[i[k]]+'.'+suffix+'.psd'
        #Compute spectra
        fn, npsd, nu = tsa.multi_taper_psd(n[0].data,Fs=1./n[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
        fe, epsd, nu = tsa.multi_taper_psd(e[0].data,Fs=1./e[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
        fu, upsd, nu = tsa.multi_taper_psd(u[0].data,Fs=1./u[0].stats.delta,adaptive=True,jackknife=False,low_bias=True)
        #Convert to dB
        npsd=10*log10(npsd)
        epsd=10*log10(epsd)
        upsd=10*log10(upsd)
        #Write to file
        savez(outpath+outname,fn=fn,fe=fe,fu=fu,npsd=npsd,epsd=epsd,upsd=upsd)



def fakequakes_allpsd(home,project_name,run_name):
    '''
    Compute PSDs for either all the synthetics fo a aprticular run or all the data
    '''
    
    from numpy import savez
    from obspy import read
    import nitime.algorithms as tsa
    from glob import glob
    from os import path,makedirs
    
    #Decide what I'm going to work on
    paths=glob(home+run_name+'/output/waveforms/'+run_name+'.*')
    for k in range(len(paths)):
        waveforms=glob(paths[k]+'/*.sac')
        print('Working on '+paths[k])
        outpath=home+project_name+'/analysis/frequency/'+paths[k].split('/')[-1]
        if not path.exists(outpath):
            makedirs(outpath)
        for ksta in range(len(waveforms)):
            sta=waveforms[ksta].split('/')[-1].split('.')[0]
            n=read(paths[k]+'/'+sta+'.LYN.sac')
            e=read(paths[k]+'/'+sta+'.LYE.sac')
            u=read(paths[k]+'/'+sta+'.LYZ.sac')
            outname=sta+'.psd'
            #Compute spectra
            fn, npsd, nu = tsa.multi_taper_psd(n[0].data,Fs=1./n[0].stats.delta,adaptive=True,jackknife=False,low_bias=True,NFFT=512)
            fe, epsd, nu = tsa.multi_taper_psd(e[0].data,Fs=1./e[0].stats.delta,adaptive=True,jackknife=False,low_bias=True,NFFT=512)
            fu, upsd, nu = tsa.multi_taper_psd(u[0].data,Fs=1./u[0].stats.delta,adaptive=True,jackknife=False,low_bias=True,NFFT=512)
            #Write to file
            
            savez(outpath+'/'+outname,fn=fn,fe=fe,fu=fu,npsd=npsd,epsd=epsd,upsd=upsd)



            
def allcoherence(home,project_name,run_name,run_number,GF_list,v_or_d,decimate,lowpass):
    '''
    Compute PSDs for either all the synthetics fo a aprticular run or all the data
    '''
    
    from numpy import genfromtxt,where,savez,c_,pi
    from obspy import read
    from mudpy.forward import lowpass as lfilter
    from mudpy.green import stdecimate 
    import nitime.algorithms as tsa
    
    #Regularized coherence
    eps=0.000001
    al=10.
    #Decide what I'm going to work on
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[4,5],dtype='f')
    datapath=home+project_name+'/data/waveforms/'
    synthpath=home+project_name+'/output/inverse_models/waveforms/'
    outpath=home+project_name+'/analysis/frequency/'
    if v_or_d.lower()=='d':
        kgf=0 #disp
        datasuffix='kdisp'
        synthsuffix='disp'
    elif v_or_d.lower()=='v':
        kgf=1 #disp
        datasuffix='kvel'
        synthsuffix='vel'
    i=where(gf[:,kgf]==1)[0]
    for k in range(len(i)):
        print('Working on '+sta[i[k]])
        #Read data
        n=read(datapath+sta[i[k]]+'.'+datasuffix+'.n')
        e=read(datapath+sta[i[k]]+'.'+datasuffix+'.e')
        u=read(datapath+sta[i[k]]+'.'+datasuffix+'.u')
        if lowpass!=None:
            fsample=1./e[0].stats.delta
            e[0].data=lfilter(e[0].data,lowpass,fsample,10)
            n[0].data=lfilter(n[0].data,lowpass,fsample,10)
            u[0].data=lfilter(u[0].data,lowpass,fsample,10)
        if decimate!=None:
            n[0]=stdecimate(n[0],decimate)
            e[0]=stdecimate(e[0],decimate)
            u[0]=stdecimate(u[0],decimate)
        #Read synthetics
        nsyn=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.n.sac')
        esyn=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.e.sac')
        usyn=read(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.u.sac')
        #What's the sampling rate correction?
        Fs=1./n[0].stats.delta
        Fsc=Fs/(2.*pi)
        #Compute coherence
        data=c_[n[0].data,nsyn[0].data].T
        fn,cn=tsa.cohere.coherence_regularized(data,epsilon=eps,alpha=al,
            csd_method={'this_method':'multi_taper_csd','adaptive':True,'low_bias':True})
        fn=fn*Fsc
        cn=cn[1,0,:].real
        data=c_[e[0].data,esyn[0].data].T
        fe,ce=tsa.cohere.coherence_regularized(data,epsilon=eps,alpha=al,
            csd_method={'this_method':'multi_taper_csd','adaptive':True,'low_bias':True})
        fe=fe*Fsc
        ce=ce[1,0,:].real
        data=c_[u[0].data,usyn[0].data].T
        fu,cu=tsa.cohere.coherence_regularized(data,epsilon=eps,alpha=al,
            csd_method={'this_method':'multi_taper_csd','adaptive':True,'low_bias':True})
        fu=fu*Fsc
        cu=cu[1,0,:].real
        #Write to file
        outname=run_name+'.'+run_number+'.'+sta[i[k]]+'.'+synthsuffix+'.coh'
        savez(outpath+outname,fn=fn,fe=fe,fu=fu,cn=cn,ce=ce,cu=cu) 
        
def source_spectra(home,project_name,run_name,run_number,rupt,nstrike,ndip):
    '''
    Tile plot of subfault source-time functions
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,savez,mean
    from mudpy.forward import get_source_time_function,add2stf
    import nitime.algorithms as tsa
    
    outpath=home+project_name+'/analysis/frequency/'
    f=genfromtxt(rupt)
    num=f[:,0]
    nfault=nstrike*ndip
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
    all=zeros(len(all_ss)*2)
    iss=2*arange(0,len(all)/2,1)
    ids=2*arange(0,len(all)/2,1)+1
    all[iss]=all_ss
    all[ids]=all_ds
    #Now parse for multiple rupture speeds
    unum=unique(num)
    #Count number of windows
    nwin=len(where(num==unum[0])[0])
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    for kfault in range(nfault):
        if kfault%10==0:
            print('... working on subfault '+str(kfault)+' of '+str(nfault))
        #Get rupture times for subfault windows
        i=where(num==unum[kfault])[0]
        trup=f[i,12]
        #Get slips on windows
        ss=all_ss[i]
        ds=all_ds[i]
        #Add it up
        slip=(ss**2+ds**2)**0.5
        #Get first source time function
        t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
        #Convert to slip rate
        s=M1/(mu[kfault]*area[kfault])
        #remove mean
        s=s-mean(s)
        #Done now compute spectra of STF
        fsample=1./(t1[1]-t1[0])
        freq, psd, nu = tsa.multi_taper_psd(s,Fs=fsample,adaptive=True,jackknife=False,low_bias=True)
        outname=run_name+'.'+run_number+'.sub'+str(kfault).rjust(4,'0')+'.stfpsd'
        savez(outpath+outname,freq=freq,psd=psd) 
        