def subfault_STFs(rupt,epicenter,nstrike,ndip,beta=None,covfile=None):
    '''
    Extract subfault source-time functions
    
    If analyzing an output .inv file make beta=0
    '''
    from numpy import genfromtxt,unique,zeros,where,meshgrid,linspace,load,arange,expand_dims,squeeze,tile,r_
    from mudpy.forward import get_source_time_function,add2stf
    from mudpy.inverse import d2epi,ds2rot
    

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
    #Compute CI
    #Load covariances
    if covfile!=None:
        rot=ds2rot(expand_dims(all,1),beta)
        C=load(covfile)
        CIplus=squeeze(rot)+1*(C**0.5)
        CIminus=squeeze(rot)-1*(C**0.5)
        CIminus[CIminus<0]=0
        slipCIplus=(CIplus[iss]**2+CIplus[ids]**2)**0.5
        slipCIminus=(CIminus[iss]**2+CIminus[ids]**2)**0.5
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
    #Get coordinates and compute distances
    source=f[0:len(unum),1:4]
    d=d2epi(epicenter,source)
    #Loop over subfaults
    Mmax=0
    Mout=[]
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
        if covfile !=None:
            slip_plus=slipCIplus[i]
            slip_minus=slipCIminus[i]
        #Get first source time function
        t1,M1=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip[0])
        if covfile !=None:
            t1plus,M1plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_plus[0])
            t1minus,M1minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[0],slip_minus[0])
        #Loop over windows
        for kwin in range(nwin-1):
            #Get next source time function
            t2,M2=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip[kwin+1])
            if covfile !=None:
                t2plus,M2plus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_plus[kwin+1])
                t2minus,M2minus=get_source_time_function(mu[kfault],area[kfault],rise_time[kfault],trup[kwin+1],slip_minus[kwin+1])
            #Add the soruce time functions
            t1,M1=add2stf(t1,M1,t2,M2)
            if covfile !=None:
                t1plus,M1plus=add2stf(t1plus,M1plus,t2plus,M2plus)
                t1minus,M1minus=add2stf(t1minus,M1minus,t2minus,M2minus)
        #Save M1 for output
        if kfault==0:
            Mout=expand_dims(M1,1).T
            tout=expand_dims(t1,1).T
        else:
            Mout=r_[Mout,expand_dims(M1,1).T]
            tout=r_[tout,expand_dims(t1,1).T]
        #Track maximum moment
        Mmax=max(Mmax,M1.max())
    print 'Maximum moment was '+str(Mmax)+'N-m'
    
    return tout,Mout
    
def fault_scaling(Mw,mu):
    '''
    Use scaling relationships of Blaser et al. 2010 to determine fault length,
    width and mean slip
    '''
    from numpy import log10
    W=10**(-1.86+0.46*Mw)
    L=10**(-2.37+0.57*Mw)
    #Mw= (2/3)log10(m0)-6.07
    M0=10**(1.5*(Mw+6.07))
    d=M0/(mu*L*1000*W*1000)
    return d,L,W
    
def make_scaling_fault(home,project_name,disp,length,width,strike,dip,hypocenter):
    '''
    Make planar fault geometry from information about fault scaling, hypocenter
    coordinates and assumed strike and dip
    '''
    
    
    
    
def make_planar_geometry(strike,dip,nstrike,dx_dip,dx_strike,epicenter,num_updip,num_downdip,rise_time):
    '''
    Make a planar fault
    
    strike - Strike angle (degs)
    dip - Dip angle (degs)
    '''
    from numpy import arange,sin,cos,deg2rad,r_,ones,arctan,rad2deg,zeros,isnan,unique,where,argsort,c_
    import pyproj
    
    proj_angle=180-strike #Angle to use for sin.cos projection (comes from strike)
    y=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    x=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    z=ones(x.shape)*epicenter[2]
    y=y*cos(deg2rad(strike))
    x=x*sin(deg2rad(strike))   
    #Save teh zero line
    y0=y.copy()
    x0=x.copy()
    z0=z.copy()
    #Initlaize temp for projection up/down dip
    xtemp=x0.copy()
    ytemp=y0.copy()
    ztemp=z0.copy()
    #Get delta h and delta z for up/ddx_dip=1own dip projection
    dh=dx_dip*cos(deg2rad(dip))
    dz=dx_dip*sin(deg2rad(dip))
    #Project updip lines
    for k in range(num_updip):
        xtemp=xtemp+dh*cos(deg2rad(proj_angle))
        ytemp=ytemp+dh*sin(deg2rad(proj_angle))
        ztemp=ztemp-dz
        x=r_[x,xtemp]
        y=r_[y,ytemp]
        z=r_[z,ztemp]
    #Now downdip lines
    xtemp=x0.copy()
    ytemp=y0.copy()
    ztemp=z0.copy()
    for k in range(num_downdip):
        xtemp=xtemp-dh*cos(deg2rad(proj_angle))
        ytemp=ytemp-dh*sin(deg2rad(proj_angle))
        ztemp=ztemp+dz
        x=r_[x,xtemp]
        y=r_[y,ytemp]
        z=r_[z,ztemp]
    #Now use pyproj to dead reckon anf get lat/lon coordinates of subfaults
    g = pyproj.Geod(ellps='WGS84')
    #first get azimuths of all points, go by quadrant
    az=zeros(x.shape)
    for k in range(len(x)):
        if x[k]>0 and y[k]>0:
            az[k]=rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]>0:
            az[k]=360+rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
        if x[k]>0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
    #Quadrant correction
    #Now horizontal distances
    d=((x**2+y**2)**0.5)*1000
    #Now reckon
    lo=zeros(len(d))
    la=zeros(len(d))
    for k in range(len(d)):
        if isnan(az[k]): #No azimuth because I'm on the epicenter
            print 'Point on epicenter'
            lo[k]=epicenter[0]
            la[k]=epicenter[1]
        else:
            lo[k],la[k],ba=g.fwd(epicenter[0],epicenter[1],az[k],d[k]) 
    #Sort them from top right to left along dip
    zunique=unique(z)
    for k in range(len(zunique)):
        i=where(z==zunique[k])[0] #This finds all faults at a certain depth
        isort=argsort(la[i]) #This sorths them south to north
        if k==0: #First loop
            laout=la[i][isort]
            loout=lo[i][isort]
            zout=z[i][isort]
        else:
            laout=r_[laout,la[i][isort]]
            loout=r_[loout,lo[i][isort]]
            zout=r_[zout,z[i][isort]]
    #Write to file
    strike=ones(loout.shape)*strike
    dip=ones(loout.shape)*dip
    tw=ones(loout.shape)*0.5
    rise=ones(loout.shape)*rise_time
    L=ones(loout.shape)*dx_strike*1000
    W=ones(loout.shape)*dx_dip*1000
    # Make output
    fault=c_[fault_num,loout,laout,zout,strike,dip,tw,rise,L,W]
    return fault