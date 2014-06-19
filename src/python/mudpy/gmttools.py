#def make_slip_slice(rupt,out):
#    '''
#    Make .slip files for GMT animation
#    '''
#    
#    from numpy import genfromtxt,savetxt,unique,where,zeros,arange,intersect1d,trapz,c_
#    from string import rjust
#    
#    #Run parameters
#    
#    dt=15
#    cumul=0 #Time slices or cumulative slip?
#    #########
#    
#    delta_t=0.01
#    f=genfromtxt(rupt)
#    trupt=f[:,12]
#    trise=f[:,7]
#    all_ss=f[:,8]
#    all_ds=f[:,9]
#    num=f[:,0]
#    #Get other parameters
#    #lon=f[0:len(unum),1]
#    #lat=f[0:len(unum),2]
#    #strike=f[0:len(unum),4]
#    #Decide on time vector
#    tslice=arange(0,trupt.max()+dt,dt)
#    #Determine triangle height at all subfaults
#    hss=2*all_ss/trise
#    hds=2*all_ds/trise
#    #Cumulative
#    ss_cumul=zeros(len(f))
#    ds_cumul=zeros(len(f))
#    #Determine time series fo each triangle
#    t=arange(0,trupt.max()+trise[0],delta_t)
#    for kslice in range(len(tslice-2)):
#        print str(kslice)+'/'+str(len(tslice)-1)
#        #And initalize slice vectors
#        ss_slice=zeros(len(f))
#        ds_slice=zeros(len(f))
#        for kfault in range(len(f)):
#            yss=zeros(t.shape)
#            yds=zeros(t.shape)
#            #Up going
#            i1=where(t>=trupt[kfault])[0]
#            i2=where(t<=(trupt[kfault]+trise[0]/2))[0] #Ascending triangle
#            i=intersect1d(i1,i2)
#            yss[i]=(2*hss[kfault]/trise[0])*t[i]-(2*hss[kfault]*trupt[kfault]/trise[0])
#            yds[i]=(2*hds[kfault]/trise[0])*t[i]-(2*hds[kfault]*trupt[kfault]/trise[0])
#            #Down going
#            i1=where(t>(trupt[kfault]+trise[0]/2))[0]
#            i2=where(t<=(trupt[kfault]+trise[0]))[0] #Ascending triangle
#            i=intersect1d(i1,i2)
#            yss[i]=(-2*hss[kfault]/trise[0])*t[i]+(2*hss[kfault]/trise[0])*(trupt[kfault]+trise[0])
#            yds[i]=(-2*hds[kfault]/trise[0])*t[i]+(2*hds[kfault]/trise[0])*(trupt[kfault]+trise[0])
#            #Now integrate slip at pertinent time interval
#            i1=where(t>=tslice[kslice])[0]
#            i2=where(t<=tslice[kslice+1])[0]
#            i=intersect1d(i1,i2)
#            ss_slice[kfault]=trapz(yss[i],t[i])
#            ds_slice[kfault]=trapz(yds[i],t[i])
#        #Combine into single model for that time slice
#        ss_cumul=ss_cumul+ss_slice
#        ds_cumul=ds_cumul+ds_slice
#        unum=unique(num)
#        lon=f[0:len(unum),1]
#        lat=f[0:len(unum),2]
#        depth=f[0:len(unum),3]
#        ss=zeros(len(unum))
#        ds=zeros(len(unum))
#        for k in range(len(unum)):
#            if cumul==0:
#                i=where(unum[k]==num)
#                ss[k]=ss_slice[i].sum()
#                ds[k]=ds_slice[i].sum()  
#                outname='slice'  
#            else:
#                i=where(unum[k]==num)
#                ss[k]=ss_cumul[i].sum()
#                ds[k]=ds_cumul[i].sum() 
#                outname='cumul'       
#        #Write outfile
#        fname=out+rjust(str(kslice),4,'0')+'.'+outname+'.slip'
#        savetxt(fname, c_[lon,lat,depth,ss,ds],fmt='%.6f\t%.6f\t%.4f\t%.2f\t%.2f')
      
def make_total_model(rupt):
    from numpy import genfromtxt,unique,where,zeros,c_,savetxt
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    fname=rupt+'.total.slip'
    savetxt(fname, c_[lon,lat,ss,ds],fmt='%.6f\t%.6f\t%6.2f\t%6.2f')
    
    
def make_sliprate_slice(rupt,nstrike,ndip,epicenter,out,tmax,dt):
    '''
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt
    from mudpy.forward import get_source_time_function,add2stf
    from string import rjust
    
    f=genfromtxt(rupt)
    num=f[:,0]
    nfault=nstrike*ndip
    unum=unique(num)
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    depth=f[0:len(unum),3]
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
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
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    t=arange(0,tmax,dt)
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j-1
             idip[k]=i
             k+=1  
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
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
        if kfault==0: #Intialize
            M=zeros((len(t),nfault))
            T=zeros((len(t),nfault))
        Mout=interp(t,t1,M1)
        M[:,kfault]=Mout
        T[:,kfault]=t
    #Now look through time slices
    maxsr=0
    print 'Writing files...'
    for ktime in range(len(t)):
        sliprate=zeros(lon.shape)
        for kfault in range(nfault):
            i=where(T[:,kfault]==t[ktime])[0]
            sliprate[kfault]=M[i,kfault]/(mu[kfault]*area[kfault])
        maxsr=max(maxsr,sliprate.max())
        fname=out+rjust(str(ktime),4,'0')+'.sliprate'
        savetxt(fname, c_[lon,lat,depth,sliprate],fmt='%.6f\t%.6f\t%.4f\t%.6f')
    print 'Maximum slip rate was '+str(maxsr)+'m/s'
    
def make_slip_slice(rupt,nstrike,ndip,epicenter,out,tmax,dt):
    '''
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt,intersect1d
    from mudpy.forward import get_source_time_function,add2stf
    from string import rjust
    from scipy.integrate import trapz
    
    f=genfromtxt(rupt)
    num=f[:,0]
    nfault=nstrike*ndip
    unum=unique(num)
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    depth=f[0:len(unum),3]
    #Get slips
    all_ss=f[:,8]
    all_ds=f[:,9]
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
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    t=arange(0,tmax,0.1)
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j-1
             idip[k]=i
             k+=1  
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print '... working on subfault '+str(kfault)+' of '+str(nfault)
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
        if kfault==0: #Intialize
            M=zeros((len(t),nfault))
            T=zeros((len(t),nfault))
        Mout=interp(t,t1,M1)
        M[:,kfault]=Mout
        T[:,kfault]=t
    #Now look through time slices
    maxsr=0
    #Now integrate slip
    t=arange(0,tmax,dt)
    print 'Writing files...'
    for ktime in range(len(t)-1):
        slip=zeros(lon.shape)
        for kfault in range(nfault):
            i1=where(T[:,kfault]>=t[ktime])[0]
            i2=where(T[:,kfault]<t[ktime+1])[0]
            i=intersect1d(i1,i2)
            slip[kfault]=trapz(M[i,kfault]/(mu[kfault]*area[kfault]),T[i,kfault])
        maxsr=max(maxsr,slip.max())
        fname=out+rjust(str(ktime),4,'0')+'.slip'
        savetxt(fname, c_[lon,lat,depth,slip],fmt='%.6f\t%.6f\t%.4f\t%.6f')
    print 'Maximum slip rate was '+str(maxsr)+'m/s'      

