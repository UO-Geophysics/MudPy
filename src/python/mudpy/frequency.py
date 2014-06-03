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
        print 'Working on '+sta[i[k]]
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
                n=stdecimate(n,decimate)
                e=stdecimate(e,decimate)
                u=stdecimate(u,decimate)
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
    eps=0.00001
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
        print 'Working on '+sta[i[k]]
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
            n=stdecimate(n,decimate)
            e=stdecimate(e,decimate)
            u=stdecimate(u,decimate)
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