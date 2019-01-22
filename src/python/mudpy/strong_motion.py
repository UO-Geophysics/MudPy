"""
Created on Wed Dec 19 13:25:35 2018

@author: dmelgarm
"""




def envelope(n,e,z,fcorner=None,Ncomponents=3):
    '''
    Build the envelope of a 3- or 2-component accelerogram using the Hilbert
    transform as implemented in obspy.signal.filter
    '''
    
    from obspy.signal.filter import envelope
    from mudpy.forward import lowpass
    
    #remove pre-event baseline    
    n[0].data-=n[0].data[0]
    e[0].data-=e[0].data[0]
    z[0].data-=z[0].data[0]
    
    
    #Initalize per-component envelopes
    nenv=n.copy()
    eenv=e.copy()
    zenv=z.copy()
    
    #make envelopes
    nenv[0].data=envelope(n[0].data)
    eenv[0].data=envelope(e[0].data)
    zenv[0].data=envelope(z[0].data)
    
    #combine envelopes into one     
    aenv=n.copy()
    aenvf=n.copy()
    
    #How many components
    if Ncomponents==3:
        aenv[0].data=(nenv[0].data**2+eenv[0].data**2+zenv[0].data**2)**0.5
    else:
        aenv[0].data=(nenv[0].data**2+eenv[0].data**2)**0.5
    
    #Low pass filter envelope
    if fcorner==None:
        aenvf=aenv.copy()
    else:
        aenvf[0].data=lowpass(aenv[0].data,fcorner,1./aenv[0].stats.delta,2)
        
    return aenvf



def peak_envelope(n,e,z,time_interval=1.0,Ncomponents=3,method='cubic',\
                  smooth=0,univariate=False,order=3):
    '''
    Build the peak envelope of a 3- or 2-component accelerogram using spline 
    interpolation between the peaks of the total acceleration in time_interval 
    seconds
    '''
    
    from scipy.interpolate import interp1d,UnivariateSpline
    from numpy import argmax,array,r_
    
    #remove pre-event baseline    
    n[0].data-=n[0].data[0]
    e[0].data-=e[0].data[0]
    z[0].data-=z[0].data[0]
    
    
    #Initalize total accel
    a=n.copy()
    
    #How many components
    if Ncomponents==3:
        a[0].data=(n[0].data**2+e[0].data**2+z[0].data**2)**0.5
    else:
        a[0].data=(n[0].data**2+e[0].data**2)**0.5
    
    #get times and values of peaks every time_interval seconds
    Nsamples=int(time_interval/a[0].stats.delta)
    
    #loop over samples
    kstart=0
    kend=Nsamples
    k=0
    while kend<a[0].stats.npts:
        
        #ake an Nsamples slice of the times and peaks
        t=a[0].times()[kstart:kend]
        p=a[0].data[kstart:kend]
        
        #Where is the peak?
        i=argmax(p)
        
        if k==0:
            time_peak=array([t[i]])
            peaks=array([p[i]])
        else:
            time_peak=r_[time_peak,array([t[i]])]
            peaks=r_[peaks,array([p[i]])]
        
        k+=1
        kstart=kend
        kend+=Nsamples
            
    #Cubic spline interpolation
    env=a.copy()
    
    if univariate==False:
        f=interp1d(time_peak,peaks,kind=method,bounds_error=False,fill_value=0)
    else:
        
        f=UnivariateSpline(time_peak, peaks, k=order, s=smooth, ext='zeros')
    env[0].data=f(env[0].times())
        
        
    return env



def Nstack(st,N=1,normalize=True):
    '''
    Make the Nth root stack, N=1 is the simple linear stack. See McFadden et al.
    in Geophysics, 1986, equation 6
    '''

    from numpy import max,abs
    
    for k in range(len(st)):   
        
        #Normalize the data?
        if normalize==True:
            data=abs(st[k].data)/max(abs(st[k].data))
        else:
            data=abs(st[k].data)
        
        #form the stack
        if k==0:
            stack=st[0].copy()
            stack.data=data**(1./N)
        else:
            stack.data+=data**(1./N)
            
    #Finish the stack
    stack.data/=len(st)
    stack.data=stack.data**N
            
    return stack
