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
    Use scaling relationships of Blaser et al. 2010 (BSSA) to determine fault length,
    width and mean slip
    '''
    from numpy import log10
    W=10**(-1.86+0.46*Mw)
    L=10**(-2.37+0.57*Mw)
    #Mw= (2/3)log10(m0)-6.07
    M0=10**(1.5*(Mw+6.07))
    d=M0/(mu*L*1000*W*1000)
    return d,L,W
    
def make_scaling_fault(home,project_name,slip,length,width,strike,dip,rake,hypocenter,faultout,ruptout):
    '''
    Make planar fault geometry from information about fault scaling, hypocenter
    coordinates and assumed strike and dip
    '''
    from numpy import deg2rad,sin,savetxt,zeros,ones,sin,cos,deg2rad
    from mudpy.forward import get_mu
    #decide on subfault size
    nstrike=20
    ndip=10
    dx_dip=width/ndip
    dx_strike=length/nstrike
    rise_time=1
    num_updip=ndip/2
    num_downdip=ndip-num_updip-1
    #Get fault
    fault=make_planar_geometry(strike,dip,nstrike,dx_dip,dx_strike,hypocenter,num_updip,num_downdip,rise_time)
    # Figure out if fault goes too shallow
    delta_dip=dx_dip*sin(deg2rad(dip))
    minz=fault[:,3].min()
    too_shallow=True
    while too_shallow:
        if minz-delta_dip<0: #Fault breaches the surface
            print 'Fault too shallow by '+str(delta_dip-minz)+'km, adjusting down-dip width...'
            num_updip-=1
            num_downdip+=1
            #Get fault
            fault=make_planar_geometry(strike,dip,nstrike,dx_dip,dx_strike,hypocenter,num_updip,num_downdip,rise_time)
            minz=fault[:,3].min()
        else: #Fault is buried, stop iterating
            too_shallow=False
    out=zeros((len(fault),14))
    out[:,0:8]=fault[:,0:8]
    out[:,8]=slip*ones(len(fault))*cos(deg2rad(rake))
    out[:,9]=slip*ones(len(fault))*sin(deg2rad(rake))
    out[:,10:12]=fault[:,8:10]
    out[:,12]=zeros(len(fault))
    out[:,13]=30e9*ones(len(fault))
    savetxt(faultout,fault,fmt='%i\t%.6f\t%.6f\t%.3f\t%i\t%i\t%.1f\t%.1f\t%.2f\t%.2f',header='No,lon,lat,z(km),strike,dip,rise,dura,ss_len(m),ds_len(m)')
    savetxt(ruptout,out,fmt='%i\t%.6f\t%.6f\t%.3f\t%i\t%i\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
    
    
    
    
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
    fault_num=arange(len(loout))+1
    fault=c_[fault_num,loout,laout,zout,strike,dip,tw,rise,L,W]
    return fault
    
class MT:
    """
    A moment tensor class
    """
    def __init__(self, m11,m22,m33,m12,m13,m23,lon,lat,depth):
        from numpy import array
        self.tensor = array([[m11,m12,m13],[m12,m22,m23],[m13,m23,m33]])
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.mt_style='rtp'
    def moment(self):
        from numpy import sqrt
        from scipy.linalg import norm
        return norm(self.tensor/sqrt(2))
    def Mw(self):
        from numpy import log10
        Mw=(2./3)*log10(self.moment())-6.07
        return Mw
    def flip(self):
        from numpy import zeros
        if self.mt_style == 'rtp':
            M=zeros((3,3))
            M[0,0]=self.tensor[1,1]
            M[1,1]=self.tensor[2,2]
            M[2,2]=self.tensor[0,0]
            M[0,1]=-self.tensor[1,2]
            M[0,2]=self.tensor[0,1]
            M[1,2]=-self.tensor[0,2]
            M[1,0]=M[0,1]
            M[2,0]=M[0,2]
            M[2,1]=M[1,2]
            self.tensor=M
            self.mt_style='xyz'
        else:
            M=zeros((3,3))
            M[0,0]=self.tensor[2,2]
            M[1,1]=self.tensor[0,0]
            M[2,2]=self.tensor[1,1]
            M[0,1]=self.tensor[0,2]
            M[0,2]=-self.tensor[1,2]
            M[1,2]=-self.tensor[0,1]
            M[1,0]=M[0,1]
            M[2,0]=M[0,2]
            M[2,1]=M[1,2]
            self.tensor=M
            self.mt_style='rtp'  
    def eigen(self):
        from numpy.linalg import eig
        self.eig_val,self.eig_vec=eig(self.tensor)
    def compute_DC(self):
        #Bst fitting double couple
        from numpy import zeros,diag
        if self.mt_style is 'rtp':
            self.flip()
        self.eigen()
        eig1=self.eig_val[0]
        eig3=self.eig_val[2]
        self.tensor_DC=zeros((3,3))
        self.tensor_DC[0,0]=0.5*(eig1-eig3)
        self.tensor_DC[2,2]=-0.5*(eig1-eig3)
        #Rotate back to original coordinates
        #Mdct=V*Mdct*V';
        #Mdc(:,:,k)=Mdct;
        self.tensor_DC=self.eig_vec.dot(self.tensor_DC).dot(self.eig_vec.transpose())
        d=self.eig_val
        E=sum(self.eig_val)/3
        d=d-E
        self.epsilon=min(abs(d))/max(abs(d))
    def get_nodal_planes(self):
        from numpy import array,arccos,arctan2,cos,sin,pi,rad2deg,real
        from numpy.linalg import norm,eig
        self.compute_DC()
        M=self.tensor_DC
        #Get eigen vectors and values
        eig_val,eig_vec=eig(M)
        #Define new basis: tension, pressure and null axes
        if eig_val[0]>0:
            t=eig_vec[:,0]
            p=eig_vec[:,2]
        else:
            t=eig_vec[:,2]
            p=eig_vec[:,0]
        #Make sure rake of T and P axes is positive (stereonet representation)
        if t[2]<0:
            t=-t
        if p[2]<0:
            p=-p
        #Get fault normal and rake vectors
        n=0.5*(p+t)
        d=0.5*(p-t)
        #normalize
        n=n/norm(n)
        d=d/norm(d)
        #Define vectors for both fault planes, Normal vector should be pointing
        #down (positive z)
        #NP1
        if n[2]<0:
            n1=-n
            d1=-d
        else:
            n1=n
            d1=d
        #NP2
        if d[2]<0:
            n2=-d
            d2=-n
        else:
            n2=d
            d2=n
        #Compute angles
        di1=arccos(n1[2])
        st1=arctan2(n1[0],-n1[1])
        phi1=array([cos(st1),sin(st1),0])
        ra1=arccos(phi1.dot(d1))
        di2=arccos(n2[2])
        st2=arctan2(n2[0],-n2[1])
        phi2=array([cos(st2),sin(st2),0])
        ra2=arccos(phi2.dot(d2))
        #Format according to right-pinule strike and dip convention
        if st1 < 0 and st1 >= -pi:
            st1=2*pi+st1
        if st2 < 0 and st2 >= -pi:
            st2=2*pi+st2
        if st1 >= 2*pi:
            st1=st1-2*pi
        if st2 >= 2*pi:
            st2=st2-2*pi
        if ra1 < 0 and ra1 >= -pi:
            ra1=2*pi+ra1
        if ra2 < 0 and ra2 >= -pi:
            ra2=2*pi+ra2
        #Fix ambiguity in rake angle
        if d1[2]>0:  #normal faulting component present
            ra1=2*pi-ra1
        if d2[2]>0:
            ra2=2*pi-ra2
        #More formatting
        if ra1 >= 2*pi:
            ra1=ra1-2*pi
        if ra2 >= 2*pi:
            ra2=ra2-2*pi
        #Final output
        st1=rad2deg(real(st1))
        di1=rad2deg(real(di1))
        ra1=rad2deg(real(ra1))
        st2=rad2deg(real(st2))
        di2=rad2deg(real(di2))
        ra2=rad2deg(real(ra2))
        self.nodal_plane1=array([st1,di1,ra1])
        self.nodal_plane2=array([st2,di2,ra2])
    def plot(self):
        from obspy.imaging.beachball import Beachball
        Beachball([self.tensor[0,0],self.tensor[1,1],self.tensor[2,2],self.tensor[0,1],self.tensor[0,2],self.tensor[1,2]])

        
        