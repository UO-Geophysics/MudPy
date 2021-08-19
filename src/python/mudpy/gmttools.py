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
      
def make_total_model(rupt,thresh):
    from numpy import genfromtxt,unique,where,zeros,c_,savetxt,sqrt,arctan2,rad2deg
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros((len(unum),1))
    ds=zeros((len(unum),1))
    lon=f[0:len(unum),1]
    lat=f[0:len(unum),2]
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Apply threshold
    j=where(sqrt(ss**2+ds**2)<thresh)[0]
    ss[j]=0
    ds[j]=0
    fname=rupt+'.total'
    savetxt(fname, c_[f[0:len(unum),0:8],ss,ds,f[0:len(unum),10:12],f[0:len(unum),13]],fmt='%d\t%10.4f\t%10.4f\t%8.4f\t%8.2f\t%6.2f\t%6.2f\t%6.2f\t%12.4e\t%12.4e\t%8.2f\t%8.2f\t%8.4e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
    
    #Now dor akes
    rakes=rad2deg(arctan2(ds,ss))
    fname=rupt+'.rakes'
    savetxt(fname, c_[f[0:len(unum),1:3],ss,ds,rakes],fmt='%d\t%10.4f\t%10.4f\t%8.4f\t%8.2f',header='lon,lat,z(km),ss-slip(m),ds-slip(m),rake(deg)')
    
    
def make_subfault(rupt):
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
    savetxt(fname, c_[lon,lat,ss,ds],fmt='%.6f\t%.6f\t%6.2f\t%6.2f',header='# No,lon,lat,z(km),strike,dip,rise,duration(s),ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rigidity(Pa)')
    
    
    
def make_sliprate_ruptfiles(rupt,nstrike,ndip,epicenter,out,tmax,dt):
    '''
    Make sliprate .rupt files for plotting several frames with tile_slip and
    animating with ffmpeg. Only write dip slip component regardless of the actual
    rake angle. We're just plotting scalar slip rate... 
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt
    from mudpy.forward import get_source_time_function,add2stf
    
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
            print( '... working on subfault '+str(kfault)+' of '+str(nfault))
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
    #First retain original rupture information
    ruptout=f[0:len(unum),:]
    ruptout[:,8]=0 #Set SS to 0
    print( 'Writing files...')
    for ktime in range(len(t)):
        sliprate=zeros(lon.shape)
        for kfault in range(nfault):
            i=where(T[:,kfault]==t[ktime])[0]
            sliprate[kfault]=M[i,kfault]/(mu[kfault]*area[kfault])
        ruptout[:,9]=sliprate #Assign slip rate to SS component    
        maxsr=max(maxsr,sliprate.max())
        fname=out+str(ktime).rjust(4,'0')+'.sliprate'
        #Write as a rupt file
        fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
        savetxt(fname,ruptout,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
    print( 'Maximum slip rate was '+str(maxsr)+'m/s')    
    
    
def make_sliprate_slice(rupt,nstrike,ndip,epicenter,out,tmax,dt):
    '''
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt
    from mudpy.forward import get_source_time_function,add2stf
    
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
    t=arange(0,tmax,dt)
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print( '... working on subfault '+str(kfault)+' of '+str(nfault))
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
    print('Writing files...')
    for ktime in range(len(t)):
        sliprate=zeros(lon.shape)
        for kfault in range(nfault):
            i=where(T[:,kfault]==t[ktime])[0]
            sliprate[kfault]=M[i,kfault]/(mu[kfault]*area[kfault])
        maxsr=max(maxsr,sliprate.max())
        fname=out+str(ktime).rjust(4,'0')+'.sliprate'
        savetxt(fname, c_[lon,lat,depth,sliprate],fmt='%.6f\t%.6f\t%.4f\t%.6f')
    print( 'Maximum slip rate was '+str(maxsr)+'m/s')
    
   
    
     
def make_fakequakes_sliprate_slice(rupt,nstrike,ndip,epicenter,outpath,run_name,tmax,dt,stf_type='dreger'):
    '''
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt
    from mudpy.forward import get_source_time_function,add2stf,build_source_time_function
    
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
    #Get rigidities
    mu=f[0:len(unum),13]
    #Get rise times
    rise_time=f[0:len(unum),7]
    #Get areas
    area=f[0:len(unum),10]*f[0:len(unum),11]
    t=arange(0,tmax,dt)
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print( '... working on subfault '+str(kfault)+' of '+str(nfault))
        
        #Get rupture times for subfault windows
        trup=f[kfault,12]
        
        #Get slips on subfault
        ss=all_ss[kfault]
        ds=all_ds[kfault]
        
        #Add it up
        slip=(ss**2+ds**2)**0.5
        
        #Get first source time function
        if rise_time[kfault]>0:
            tstf,Mdot=build_source_time_function(rise_time[kfault],0.05,rise_time[kfault]*1.5,stf_type=stf_type,time_offset=trup)
            #What is the moment at this subfault?
            moment=slip*mu[kfault]*area[kfault]
            #Scale stf by this moment
            scale_factor=moment/dt
            Mdot=Mdot*scale_factor
        else:
            tstf=arange(2)
            Mdot=zeros(2)
        if kfault==0: #Intialize
            M=zeros((len(t),nfault))
            T=zeros((len(t),nfault))
        
        #interpolate to requested output dt
        Mout=interp(t,tstf,Mdot)
        
        #place in output vector
        M[:,kfault]=Mout
        T[:,kfault]=t
        
    #Now look through time slices
    maxsr=0
    print( 'Writing files...')
    for ktime in range(len(t)):
        sliprate=zeros(lon.shape)
        for kfault in range(nfault):
            sliprate[kfault]=M[ktime,kfault]/(mu[kfault]*area[kfault])
        maxsr=max(maxsr,sliprate.max())
        fname=outpath+run_name+'.'+str(ktime).rjust(4,'0')+'.sliprate'
        print( sliprate.max())
        savetxt(fname, c_[lon,lat,depth,sliprate],fmt='%.6f\t%.6f\t%.4f\t%.6f')
    print( 'Maximum slip rate was '+str(maxsr)+'m/s' )  
       
       
        
            
def make_usgs_sliprate_slice(rupt,outpath,tmax,dt):
    '''
    '''
    from numpy import zeros,where,arange,interp,c_,savetxt,cos,pi,roll,isnan
    from mudpy.forward import build_source_time_function
    from scipy.integrate import trapz
    
    f,segment_out,segment_data_out = read_neic_param(rupt)
    nfault=len(f)

    #coordinates
    lon=f[:,1]
    lat=f[:,0]
    depth=f[:,2]

    #Get slips
    total_slip=f[:,3]/100 #to meters  
    
    #Get rise times
    tup=f[:,8]
    tdown=f[:,9]
    
    #Get rupture times for subfault windows
    trup=f[:,7]
    
    #get subfault moments
    moment=f[:,10]/1e7

    #Output time vector
    t=arange(0,tmax,dt)
    
    #Loop over subfaults
    for kfault in range(nfault):
        if kfault%10==0:
            print('... working on subfault '+str(kfault)+' of '+str(nfault))
        
        #Add it up
        slip=total_slip[kfault]
        
        #Build local sliprate function
        s1=(1./(tup[kfault]+tdown[kfault]))*(1-cos((pi*t)/tup[kfault]))
        i=where(t>tup[kfault])[0]
        s1[i]=0
        s2=(1./(tup[kfault]+tdown[kfault]))*(1+cos((pi*(t-tup[kfault]))/tdown[kfault]))
        i=where(t<=tup[kfault])[0]
        s2[i]=0 
        i=where(t>tup[kfault]+tdown[kfault])[0]
        s2[i]=0
        #add the two 
        s=s1+s2
        #shift to right onset time
        dN=int(trup[kfault]/dt)
        s=roll(s,dN)
        #rescale to the correct slip
        area=trapz(s,t)
        scale=slip/area
        s=s*scale

        if kfault==0: #Intialize outout vectors
            S=zeros((len(t),nfault))
            T=zeros((len(t),nfault))
        
        
        #place in output vector
        S[:,kfault]=s
        T[:,kfault]=t
        
    #Now look through time slices
    maxsr=0
    print('Writing files...')
    for ktime in range(len(t)):

        sliprate=S[ktime,:]
        i=where(isnan(sliprate)==True)[0]
        sliprate[i]=0
        maxsr=max(maxsr,sliprate.max())
        fname=outpath+str(ktime).rjust(4,'0')+'.sliprate'
        print(sliprate.max())
        savetxt(fname, c_[lon,lat,depth,sliprate],fmt='%.6f\t%.6f\t%.4f\t%.6f')
    print('Maximum slip rate was '+str(maxsr)+'m/s')    
    
    #also output total slip
    out=c_[lon,lat,total_slip]   
    fname=outpath+'_total_slip.txt'
    savetxt(fname, out,fmt='%.6f\t%.6f\t%.4f')     
                    
                        
                            
                                    
def make_slip_slice(rupt,nstrike,ndip,epicenter,out,tmax,dt):
    '''
    '''
    from numpy import genfromtxt,unique,zeros,where,arange,interp,c_,savetxt,intersect1d
    from mudpy.forward import get_source_time_function,add2stf
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
    print(nfault)
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
        if kfault==0: #Intialize
            M=zeros((len(t),nfault))
            T=zeros((len(t),nfault))
        Mout=interp(t,t1,M1)
        M[:,kfault]=Mout
        T[:,kfault]=t
    #Now look through time slices
    maxsr=0
    #Now integrate slip
    t=arange(0,tmax+dt*0.1,dt)
    print(t)
    #First retain original rupture information
    ruptout=f[0:len(unum),:]
    ruptout[:,8]=0 #Set SS to 0
    print('Writing files...')
    for ktime in range(len(t)-1):
        slip=zeros(lon.shape)
        for kfault in range(nfault):
            i1=where(T[:,kfault]>=t[ktime])[0]
            i2=where(T[:,kfault]<t[ktime+1])[0]
            i=intersect1d(i1,i2)
            slip[kfault]=trapz(M[i,kfault]/(mu[kfault]*area[kfault]),T[i,kfault])
        ruptout[:,9]=slip #Assign slip to SS component  
        maxsr=max(maxsr,slip.max())
        fname=out+str(ktime).rjust(4,'0')+'.slip'
        fmtout='%6i\t%.4f\t%.4f\t%8.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%12.4e\t%12.4e%10.1f\t%10.1f\t%8.4f\t%.4e'
        savetxt(fname,ruptout,fmtout,header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
    print('Maximum slip was '+str(maxsr)+'m')      




def make_psvelo(stafile,directory,fout,run_name,run_number):
    '''
    Make psvelo file
    '''
    from numpy import genfromtxt,zeros,savetxt
    from glob import glob
    
    
    sta=genfromtxt(stafile,usecols=0,dtype='S')
    out=zeros((len(sta),7))
    lonlat=genfromtxt(stafile,usecols=[1,2],dtype='S')
    for k in range(len(sta)):
        out[k,0:2]=lonlat[k,:]
        #neu=genfromtxt(glob(directory+sta[k]+'*.neu')[0])
        neu=genfromtxt(glob(directory+run_name+'.'+run_number+'.'+sta[k]+'*.neu')[0])
        neu=neu*1000
        out[k,2]=neu[1] #East
        out[k,3]=neu[0] #North
        out[k,4]=neu[2] #Up
        out[k,5]=0.001 #East sigma
        out[k,6]=0.001 #North sigma
    savetxt(fout,out,fmt='%10.4f\t%10.4f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t')
    

def final_dtopo(dtopo_file,out_file):
    '''
    Extract final displacememnt from dtopo file
    '''
    from numpy import genfromtxt,savetxt,where
    
    #Read dtopo
    dtopo=genfromtxt(dtopo_file)
    tmax=dtopo[:,0].max()
    i=where(dtopo[:,0]==tmax)[0]
    out=dtopo[i,1:]
    savetxt(out_file,out,fmt='%10.6f\t%10.6f\t%10.6f')

def insar_xyz(home,project_name,run_name,run_number,GF_list,outfile):
    '''
    Make an xyz file for plotting InSAR sub_sampled residuals
    '''
    from numpy import genfromtxt,where,zeros,c_,savetxt
    
    #Decide what to plot
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='U')
    lon_all=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1],dtype='f')
    lat_all=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[2],dtype='f')
    gf=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[7],dtype='f')
    datapath=home+project_name+'/data/statics/'
    synthpath=home+project_name+'/output/inverse_models/statics/'
    #synthpath=home+project_name+'/output/forward_models/'
    i=where(gf==1)[0] #Which stations have statics?
    lon=lon_all[i]
    lat=lat_all[i]
    los_data=zeros(len(i))
    los_synth=zeros(len(i))
    for k in range(len(i)):
        neu=genfromtxt(datapath+sta[i[k]]+'.los')
        #neu=genfromtxt(datapath+sta[i[k]]+'.static.neu')
        los_data[k]=neu[0]
        neus=genfromtxt(synthpath+run_name+'.'+run_number+'.'+sta[i[k]]+'.los')
        #neus=genfromtxt(synthpath+sta[i[k]]+'.static.neu')
        los_synth[k]=neus[0]
    #Make plot  
    out=c_[lon,lat,los_data,los_synth,los_data-los_synth]  
    savetxt(outfile,out,fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f')


def make_shakemap_slice(home,project_name,run_name,time_epi,GF_list,dt,tmax):
    '''
    Make xyz files with current ground velocity
    '''
    from numpy import genfromtxt,arange,sqrt,zeros,where,c_,savetxt
    from obspy import read
    from mudpy.forward import lowpass
    
    t=arange(0,tmax,dt)
    fcorner=0.5
    sta=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=0,dtype='S')
    lonlat=genfromtxt(home+project_name+'/data/station_info/'+GF_list,usecols=[1,2])
    for ksta in range(len(sta)):
        n=read(home+project_name+'/output/forward_models/'+run_name+'.'+sta[ksta]+'.vel.n')
        e=read(home+project_name+'/output/forward_models/'+run_name+'.'+sta[ksta]+'.vel.n')
        n[0].data=lowpass(n[0].data,fcorner,1./n[0].stats.delta,2)
        e[0].data=lowpass(n[0].data,fcorner,1./e[0].stats.delta,2)
        if ksta==0:
            h=n.copy()
        else:
            h+=n[0].copy()
        h[ksta].data=sqrt(n[0].data**2+e[0].data**2)
        h[ksta].trim(starttime=time_epi)
        print(h[ksta].stats.starttime)
    vout=zeros(len(lonlat))
    maxv=0
    for kt in range(len(t)):
        for ksta in range(len(sta)):
            i=where(h[ksta].times()==t[kt])[0]
            vout[ksta]=h[ksta].data[i]
        if vout.max()>maxv:
            maxv=vout.max()
        out=c_[lonlat,vout]
        num=str(kt).rjust(4,'0')
        print(num)
        savetxt(home+project_name+'/analysis/shake/'+num+'.shake',out,fmt='%10.6f\t%10.6f\t%10.4f')
    print('Max velocity was '+str(maxv)+'m/s')
        


    
def read_neic_param(fault_file):
    '''
    Parse the param text format
    '''
    
    
    from numpy import r_,c_,array,expand_dims
    
    f=open(fault_file)
    #get number of segments
    line=f.readline()
    
    fault_out=array([])
    segment_out=[]
    segment_data_out=[]
    kread=0
    while True:
        
        #read fault segment header
        line=f.readline()
        if line=='':
            break
        Dx=float(line.split()[6].replace('km',''))*1000
        Dy=float(line.split()[10].replace('km',''))*1000
        Nsubfaults=int(line.split()[4])*int(line.split()[8])
        
        #Read fault segment boundary information
        line=f.readline() #skip two lines
        line=f.readline()
        for k in range(5):
            line=f.readline()
            if k==0:
                segment=expand_dims(array(line.split()).astype('float'),0)
            else:
                segment=r_[segment,expand_dims(array(line.split()).astype('float'),0)]
        #Append segment boundary
        segment_out.append(segment)
        line=f.readline() #skip a line
        
        #read segment subfault data
        for k in range(Nsubfaults):
            line=f.readline()
            
            if kread==0:
                fault_out=c_[expand_dims(array(line.split()).astype('float'),0),Dx,Dy]
            else:
                fault_out=r_[fault_out,c_[expand_dims(array(line.split()).astype('float'),0),Dx,Dy]]   
            kread+=1
            #Save segment data
            if k==0:
                segment_data=expand_dims(array(line.split()).astype('float'),0)
            else:
                segment_data=r_[segment_data,expand_dims(array(line.split()).astype('float'),0)]  
        #Append
        segment_data_out.append(segment_data)
            
        #Done
        if line=='':
            break
    
    return fault_out,segment_out,segment_data_out        
                


def triangular_rupt_2_gmt(meshfile,slipfile,outfile,kinematic_out_folder=None,percentage=0):
    
    '''
    DM Note: Modified from Brendan's script because he refused to do a pull request :)
    
    This will take a triangular mesh and a mudpy slip file and create a total slip file
    hat can be plotted in GMT and time slices and single seconds to plot the kinematic
    #slip model
    Written by Brendan Crowell, July 30, 2020
    
    IN:
        meshfile = 'simeonof_slab2.mshout' #location of triangular mesh file
        slipfile = 'final_10patch.0001.inv' #slip model with subevents
        out_file = 'test.txt' #the GMT .xy file to plot with psxy
        kinslipfolder = 'kinematic/' #create this in your working folder
            
    OUT:
        Nothing
    
    '''
    
    import numpy
    import csv
    import math
    
    
    #Read mesh file first
    
    meshnum = list()
    meshlon1 = list()
    meshlon2 = list()
    meshlon3 = list()
    meshlat1 = list()
    meshlat2 = list()
    meshlat3 = list()
    meshdep1 = list()
    meshdep2 = list()
    meshdep3 = list()
    faultarea = list()
    with open(meshfile, 'r') as f:
        next(f)
        reader = csv.reader(f,delimiter='\t')
        for row in reader:
            meshnum.append(int(row[0]))
            meshlon1.append(float(row[4]))
            meshlat1.append(float(row[5]))
            meshdep1.append(float(row[6]))
            meshlon2.append(float(row[7]))
            meshlat2.append(float(row[8]))
            meshdep2.append(float(row[9]))
            meshlon3.append(float(row[10]))
            meshlat3.append(float(row[11]))
            meshdep3.append(float(row[12]))
            faultarea.append(float(row[14]))
    
    
    #Read inverse file
    
    invnum = list()
    ss = list()
    ds = list()
    rupttime = list()
    risetime = list()
    duration = list()
    rig = list()
    
    
    with open(slipfile, 'r') as f:
        next(f)
        reader = csv.reader(f,delimiter='\t')
        for row in reader:
            invnum.append(int(row[0]))
            ss.append(float(row[8]))
            ds.append(float(row[9]))
            #rupttime.append(float(row[12]))
            risetime.append(float(row[6]))
            duration.append(float(row[7]))
            rig.append(float(row[13]))
    
    
    INVN = numpy.asarray(invnum)
    MESN = numpy.asarray(meshnum)
    
    SS = numpy.asarray(ss)
    DS = numpy.asarray(ds)
    
    RT = numpy.asarray(rupttime)
    DR = numpy.asarray(duration)
    
    TOTS = numpy.sqrt(numpy.power(SS,2)+numpy.power(DS,2))
    FA = numpy.asarray(faultarea)
    RIG = numpy.asarray(rig)
    
    
    #Total slip model
    moment = 0
    fso = open(outfile,'w')
    slip_threshold=(percentage/100)*TOTS.max()
    for i in range(0, numpy.amax(MESN)):
        a1 = numpy.where(MESN[i] == INVN)[0]
        totslip = numpy.sum(TOTS[a1])
#        print (i+1,totslip*100)
        if (totslip >= slip_threshold):
            moment = moment+FA[i]*1000*1000*numpy.mean(RIG[a1])*totslip
            lon1 = "{0:.4f}".format(meshlon1[i])
            lon2 = "{0:.4f}".format(meshlon2[i])
            lon3 = "{0:.4f}".format(meshlon3[i])
            lat1 = "{0:.4f}".format(meshlat1[i])
            lat2 = "{0:.4f}".format(meshlat2[i])
            lat3 = "{0:.4f}".format(meshlat3[i])
            dep1 = "{0:.4f}".format(meshdep1[i])
            dep2 = "{0:.4f}".format(meshdep2[i])
            dep3 = "{0:.4f}".format(meshdep3[i])
            ts = "{0:.4f}".format(totslip)
            fso.write('> -Z'+ts+'\n')
            fso.write(lon1+' '+lat1+' '+dep1+'\n')
            fso.write(lon2+' '+lat2+' '+dep2+'\n')
            fso.write(lon3+' '+lat3+' '+dep3+'\n')
            fso.write(lon1+' '+lat1+' '+dep1+'\n')
    
    fso.close()
#    print(moment,(math.log10(moment)-9.1)/1.5)



def gmtColormap(fileName):
      '''
      Convert a cpt GMT color palette file into a matplotlib colormap object.
      Note that not all default GMT palettes can be converted with this method.
      For example hot.cpt will fail but seis.cpt will be fine. If you need
      more cpt files checkout cpt-city (http://soliton.vm.bytemark.co.uk/pub/cpt-city/)
      
      Usage:
          cm = gmtColormap(fileName)
        
      IN:
          fileName: Absolute path to cpt file
      OUT:
          cm: matplotlib colormap obejct
          
      Example:
          
          colorMap = gmtColormap('/opt/local/share/gmt/cpt/haxby.cpt')
          
      you can then use your colorMap object as you would any other matplotlib
      object, for example:
         
         pcolor(x,y,cmap=colorMap)
          
      
      Diego Melgar 2,2014, modified from original by James Boyle
      '''
      
      from matplotlib import colors
      import colorsys
      import numpy as N

      f = open(fileName)
      lines = f.readlines()
      f.close()

      x = []
      r = []
      g = []
      b = []
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
             pass
          else:
              x.append(float(ls[0]))
              r.append(float(ls[1]))
              g.append(float(ls[2]))
              b.append(float(ls[3]))
              xtemp = float(ls[4])
              rtemp = float(ls[5])
              gtemp = float(ls[6])
              btemp = float(ls[7])

      x.append(xtemp)
      r.append(rtemp)
      g.append(gtemp)
      b.append(btemp)

      x = N.array( x )
      r = N.array( r )
      g = N.array( g )
      b = N.array( b )
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "RGB":
          r = r/255.
          g = g/255.
          b = b/255.
      xNorm = (x - x[0])/(x[-1] - x[0])

      red = []
      blue = []
      green = []
      for i in range(len(x)):
          red.append([xNorm[i],r[i],r[i]])
          green.append([xNorm[i],g[i],g[i]])
          blue.append([xNorm[i],b[i],b[i]])
      colorDict = {"red":red, "green":green, "blue":blue}
      cmap=colors.LinearSegmentedColormap('cmap',colorDict,256)
      return (cmap)
    